import type { FileRef } from '$lib/schema/model'
import type { DiatomicsCurves, ModelData } from '$lib/types'
import { readFile } from 'node:fs/promises'
import { isAbsolute, resolve } from 'node:path'
import { setTimeout as sleep } from 'node:timers/promises'
import { gunzipSync } from 'node:zlib'

const repo_root = resolve(import.meta.dirname, `../../../..`)

type DiatomicsSource = { pred_file?: FileRef | null }
type FetchOptions = { fetch_fn?: typeof fetch; root_dir?: string; max_attempts?: number }

const parse_gzipped_json = (
  bytes: Uint8Array | ArrayBuffer,
  source: string,
): DiatomicsCurves => {
  const uint8_bytes = bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes)
  if (uint8_bytes.length === 0) throw new Error(`${source} returned empty response body`)
  if (uint8_bytes[0] !== 0x1f || uint8_bytes[1] !== 0x8b) {
    throw new Error(`${source} did not return gzip-compressed JSON`)
  }
  return JSON.parse(gunzipSync(uint8_bytes).toString(`utf-8`)) as DiatomicsCurves
}

// Missing cache files may be downloaded; other filesystem errors must propagate.
const read_local_prediction = (path: string): Promise<Uint8Array | undefined> =>
  readFile(path).catch((error: unknown) => {
    if (!(error instanceof Error && `code` in error && error.code === `ENOENT`)) {
      throw error
    }
    return undefined
  })

const prediction_path = (name: string, root_dir: string): string =>
  isAbsolute(name) ? name : resolve(root_dir, name)

// the canonical figshare.com/files/<id> pred_file.url sits behind an AWS WAF "challenge"
// that 202s any non-browser (server-side) request; rewrite it to the ndownloader host,
// which serves the identical file via a plain signed-S3 redirect with no challenge
const to_download_url = (url: string): string =>
  url.replace(
    /^https:\/\/figshare\.com\/files\//,
    `https://ndownloader.figshare.com/files/`,
  )

const fetch_remote_diatomics_once = async (
  pred_file_url: string,
  fetch_fn: typeof fetch,
): Promise<DiatomicsCurves> => {
  const controller = new AbortController()
  const timeout = setTimeout(() => controller.abort(), 30_000)
  try {
    const response = await fetch_fn(to_download_url(pred_file_url), {
      signal: controller.signal,
    })
    const waf_action = response.headers.get(`x-amzn-waf-action`)
    if (waf_action) throw new Error(`Figshare WAF challenge: ${waf_action}`)
    if (response.status !== 200) {
      throw new Error(`${response.status} ${response.statusText}`)
    }
    return parse_gzipped_json(await response.arrayBuffer(), pred_file_url)
  } finally {
    clearTimeout(timeout)
  }
}

export async function fetch_diatomics_data(
  { pred_file }: DiatomicsSource,
  { fetch_fn = fetch, root_dir = repo_root, max_attempts = 3 }: FetchOptions = {},
): Promise<DiatomicsCurves> {
  if (pred_file?.name) {
    const path = prediction_path(pred_file.name, root_dir)
    const bytes = await read_local_prediction(path)
    if (bytes) return parse_gzipped_json(bytes, path)
  }

  const pred_file_url = pred_file?.url
  if (!pred_file_url) throw new Error(`No local diatomics file or remote URL`)

  if (!Number.isInteger(max_attempts) || max_attempts < 1) {
    throw new Error(`max_attempts must be a positive integer, got ${max_attempts}`)
  }
  const fetch_with_retry = async (attempt = 0): Promise<DiatomicsCurves> => {
    try {
      return await fetch_remote_diatomics_once(pred_file_url, fetch_fn)
    } catch (error) {
      if (error instanceof Error && error.message.startsWith(`Figshare WAF challenge:`)) {
        throw error
      }
      if (attempt + 1 >= max_attempts) {
        throw new Error(
          `${pred_file_url} failed after ${max_attempts} attempts: ${String(error)}`,
          { cause: error },
        )
      }
      await sleep(1000 * 2 ** attempt)
      return fetch_with_retry(attempt + 1)
    }
  }

  return fetch_with_retry()
}

type MdPerSystemRow = Record<string, number | string>

// per-system energy/force RMSEs are stored in eV; report meV to match the model-level
// MD table units
const mev_cols = new Set([`energy_rmse`, `force_rmse`])

// Read local predictions or their published URL. Unpublished, unavailable files
// return null; malformed files and failed downloads must surface during the build.
export async function read_md_per_system(
  model: Pick<ModelData, `metrics`>,
  root_dir: string = repo_root,
): Promise<MdPerSystemRow[] | null> {
  const prediction = model.metrics?.md?.pred_file
  if (!prediction) return null

  let bytes = await read_local_prediction(prediction_path(prediction.name, root_dir))
  if (!bytes) {
    if (!prediction.url) return null
    const url = to_download_url(prediction.url)
    const response = await fetch(url, { signal: AbortSignal.timeout(30_000) })
    if (response.status !== 200) {
      throw new Error(`MD predictions ${url}: HTTP ${response.status}`)
    }
    bytes = new Uint8Array(await response.arrayBuffer())
  }

  const csv = new TextDecoder()
    .decode(bytes[0] === 0x1f && bytes[1] === 0x8b ? gunzipSync(bytes) : bytes)
    .replace(/^\uFEFF/, ``)
  const rows: string[][] = []
  let fields: string[] = []
  let consumed = 0
  // Sticky matching consumes whole fields, including quoted commas, escaped quotes,
  // and embedded newlines. A gap means malformed CSV, never a skipped field.
  for (const [token, value, delimiter] of csv.matchAll(
    /(?<value>"(?:[^"]|"")*"|[^",\r\n]*)(?<delimiter>,|\r\n|\r|\n|$)/gy,
  )) {
    if (!token && !fields.length) break
    fields.push(value.startsWith(`"`) ? value.slice(1, -1).replaceAll(`""`, `"`) : value)
    consumed += token.length
    if (delimiter !== `,`) {
      rows.push(fields)
      fields = []
    }
  }
  const [cols, ...data] = rows
  if (
    consumed !== csv.length ||
    !cols?.includes(`system`) ||
    cols.some((col) => !col) ||
    new Set(cols).size !== cols.length ||
    data.some((row) => row.length !== cols.length)
  ) {
    throw new Error(`Invalid MD prediction CSV: ${prediction.name}`)
  }
  return data.map((values) => {
    const row: MdPerSystemRow = {}
    for (const [idx, raw] of values.entries()) {
      if (raw === ``) continue // missing value (e.g. stress-less system pressure)
      const num = Number(raw)
      row[cols[idx]] = Number.isFinite(num)
        ? num * (mev_cols.has(cols[idx]) ? 1000 : 1)
        : raw
    }
    return row
  })
}
