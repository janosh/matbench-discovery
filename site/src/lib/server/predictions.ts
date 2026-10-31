import type { FileRef } from '#lib/schema/model.js'
import type { DiatomicsCurves, ModelData } from '#lib/types.js'
import { createHash } from 'node:crypto'
import { readFile } from 'node:fs/promises'
import { isAbsolute, resolve } from 'node:path'
import { setTimeout as sleep } from 'node:timers/promises'
import { gunzipSync } from 'node:zlib'

// vite and vitest run from site/. import.meta.dirname breaks once bundled: prerender runs
// this from .svelte-kit/output/server/chunks/, which made every build re-download local
// prediction files from Figshare.
const repo_root = resolve(process.cwd(), `..`)

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

// Retrying cannot fix a WAF challenge or a 4xx (missing file, no access)
class PermanentFetchError extends Error {
  constructor(message: string) {
    super(message)
    this.name = 'PermanentFetchError'
  }
}

// Download one prediction file's bytes from its published URL. Figshare intermittently
// serves 5xx or drops connections, so those retry with exponential backoff.
const fetch_prediction_bytes = async (
  url: string,
  fetch_fn: typeof fetch,
  max_attempts: number,
): Promise<Uint8Array> => {
  if (!Number.isInteger(max_attempts) || max_attempts < 1) {
    throw new Error(`max_attempts must be a positive integer, got ${max_attempts}`)
  }
  const download_url = to_download_url(url)
  for (let attempt = 1; ; attempt++) {
    try {
      const response = await fetch_fn(download_url, {
        signal: AbortSignal.timeout(30_000),
      })
      const waf_action = response.headers.get(`x-amzn-waf-action`)
      if (waf_action)
        throw new PermanentFetchError(`Figshare WAF challenge: ${waf_action}`)
      if (response.status !== 200) {
        const msg = `HTTP ${response.status} ${response.statusText}`
        throw response.status < 500
          ? new PermanentFetchError(`${download_url}: ${msg}`)
          : new Error(msg)
      }
      return new Uint8Array(await response.arrayBuffer())
    } catch (error) {
      if (error instanceof PermanentFetchError) throw error
      if (attempt >= max_attempts) {
        const msg = `${url} failed after ${max_attempts} attempts: ${String(error)}`
        throw new Error(msg, { cause: error })
      }
      await sleep(1000 * 2 ** (attempt - 1))
    }
  }
}

// size/md5 mismatch of bytes vs their registered FileRef, undefined when they match
const registered_mismatch = (bytes: Uint8Array, file: FileRef): string | undefined => {
  if (file.size != null && bytes.length !== file.size)
    return `size ${bytes.length} != registered ${file.size}`
  const md5 = file.md5 && createHash(`md5`).update(bytes).digest(`hex`)
  return md5 && md5 !== file.md5 ? `md5 ${md5} != registered ${file.md5}` : undefined
}

// Bytes of a registered prediction file: the local copy if it matches the FileRef, else
// the verified download (a stale local file must never stand in for the registered
// artifact). undefined when there is neither a local copy nor a URL.
const read_registered_prediction = async (
  file: FileRef,
  { fetch_fn = fetch, root_dir = repo_root, max_attempts = 3 }: FetchOptions,
): Promise<{ bytes: Uint8Array; source: string } | undefined> => {
  const path = prediction_path(file.name, root_dir)
  const local = await read_local_prediction(path)
  const mismatch = local && registered_mismatch(local, file)
  if (local && !mismatch) return { bytes: local, source: path }
  if (!file.url) {
    if (mismatch)
      throw new Error(
        `Local ${path} differs from its registered artifact (${mismatch}) and no url is registered`,
      )
    return undefined
  }
  const bytes = await fetch_prediction_bytes(file.url, fetch_fn, max_attempts)
  const remote_mismatch = registered_mismatch(bytes, file)
  if (remote_mismatch)
    throw new Error(
      `${file.url} differs from its registered artifact (${remote_mismatch})`,
    )
  return { bytes, source: file.url }
}

export async function fetch_diatomics_data(
  { pred_file }: DiatomicsSource,
  options: FetchOptions = {},
): Promise<DiatomicsCurves> {
  const prediction = pred_file && (await read_registered_prediction(pred_file, options))
  if (!prediction) throw new Error(`No local diatomics file or remote URL`)
  return parse_gzipped_json(prediction.bytes, prediction.source)
}

type MdPerSystemRow = Record<string, number | string>

// per-system energy/force RMSEs are stored in eV; report meV to match the model-level
// MD table units
const mev_cols = new Set([`energy_rmse`, `force_rmse`])

// Read the registered local prediction or its published URL. Unpublished, unavailable
// files return null; malformed files and failed downloads must surface during the build.
export async function read_md_per_system(
  model: Pick<ModelData, `metrics`>,
  options: FetchOptions = {},
): Promise<MdPerSystemRow[] | null> {
  const prediction = model.metrics?.md?.pred_file
  if (!prediction) return null
  const { bytes } = (await read_registered_prediction(prediction, options)) ?? {}
  if (!bytes) return null

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
