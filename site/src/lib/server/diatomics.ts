import type { DiatomicsCurves } from '$lib/types'
import { readFile } from 'node:fs/promises'
import { dirname, isAbsolute, resolve } from 'node:path'
import { fileURLToPath } from 'node:url'
import { gunzipSync } from 'node:zlib'

const repo_root = resolve(dirname(fileURLToPath(import.meta.url)), `../../../..`)
const sleep = (ms: number): Promise<void> =>
  new Promise((resolve) => setTimeout(resolve, ms))

type DiatomicsSource = { pred_file?: string | null; pred_file_url?: string }
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

const read_local_diatomics = async (
  pred_file: string | null | undefined,
  root_dir: string,
): Promise<DiatomicsCurves | null> => {
  if (!pred_file) return null

  const file_path = isAbsolute(pred_file) ? pred_file : resolve(root_dir, pred_file)
  try {
    return parse_gzipped_json(await readFile(file_path), file_path)
  } catch (error) {
    if (error instanceof Error && `code` in error && error.code === `ENOENT`) return null
    throw error
  }
}

const fetch_remote_diatomics_once = async (
  pred_file_url: string,
  fetch_fn: typeof fetch,
): Promise<DiatomicsCurves> => {
  const controller = new AbortController()
  const timeout = setTimeout(() => controller.abort(), 30_000)
  try {
    const response = await fetch_fn(pred_file_url, { signal: controller.signal })
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
  { pred_file, pred_file_url }: DiatomicsSource,
  { fetch_fn = fetch, root_dir = repo_root, max_attempts = 3 }: FetchOptions = {},
): Promise<DiatomicsCurves> {
  const local_data = await read_local_diatomics(pred_file, root_dir)
  if (local_data) return local_data

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
