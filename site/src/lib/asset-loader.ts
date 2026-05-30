// Shared loader for gzipped JSON assets (energy-parity, kappa-parity, ...).
// Decompresses via matterviz's gzip helper, parses JSON, and caches promises.
import { decompress_data } from 'matterviz/io'

const asset_cache = new Map<string, Promise<unknown>>()

// load lifecycle status shared by the parity plot components
export type LoadStatus = `idle` | `loading` | `ready` | `error`

// per-model asset files wrap their payload as { model: ... }
export interface ModelAsset<TModel> {
  model?: TModel
}
// fields shared by every parity base asset (per-material id + formula columns)
export interface ParityBase {
  material_ids: string[]
  formulas: string[]
}
// fields shared by every per-model parity asset
export interface ParityModel {
  model_key: string
  model_label: string
}
// fields shared by every parity scatter point
export interface ParityPoint {
  material_id: string
  formula: string
}

export const clear_asset_cache = () => asset_cache.clear()

// Pick a configured base URL (e.g. from a Vite env var) over the manifest default,
// trimming trailing slashes.
export const resolve_asset_base_url = (
  configured: string | undefined,
  fallback: string,
): string => (configured?.trim() || fallback).replace(/\/+$/, ``)

export const join_asset_url = (base_url: string, asset: string): string =>
  `${base_url}/${asset.replace(/^\/+/, ``)}`

async function read_asset_text(url: string): Promise<string> {
  const response = await fetch(url)
  if (!response.ok) throw new Error(`Failed to load ${url}: HTTP ${response.status}`)

  const buffer = await response.arrayBuffer()
  const bytes = new Uint8Array(buffer)
  // raw .gz (served without Content-Encoding) starts with the gzip magic bytes and needs
  // manual decompression; if the server already decompressed it, bytes are plain JSON text
  const is_gzip = bytes[0] === 0x1f && bytes[1] === 0x8b
  return is_gzip ? decompress_data(buffer, `gzip`) : new TextDecoder().decode(bytes)
}

// Validate that a decoded asset array has the row count promised by its manifest.
// `name` should carry a caller prefix, e.g. `energy parity formulas`.
export function assert_array_length(name: string, values: unknown, length: number) {
  // describe non-arrays (e.g. a stale asset whose field is missing -> undefined) so the
  // error is clear instead of "Cannot read properties of undefined (reading 'length')"
  const got = Array.isArray(values)
    ? values.length
    : values === undefined
      ? `missing field`
      : typeof values
  if (got !== length) {
    throw new Error(`Invalid ${name}: expected ${length} rows, got ${got}`)
  }
}

export function load_json_asset<T>(url: string): Promise<T> {
  let asset = asset_cache.get(url)
  if (!asset) {
    asset = read_asset_text(url)
      .then((text) => JSON.parse(text) as T)
      .catch((error) => {
        asset_cache.delete(url)
        throw error
      })
    asset_cache.set(url, asset)
  }
  return asset as Promise<T>
}
