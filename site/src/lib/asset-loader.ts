// Shared loader for gzipped JSON assets (energy-parity, kappa-parity, ...).
// Decompresses via matterviz's gzip helper, parses JSON, and caches promises.
import { decompress_data } from 'matterviz/io'

const asset_cache = new Map<string, Promise<unknown>>()

// load lifecycle status shared by the parity plot components
export type LoadStatus = `idle` | `loading` | `ready` | `error`

// per-model asset files wrap their payload as { model: ... }
interface ModelAsset<TModel> {
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
): string => {
  const trimmed = configured?.trim() ?? ``
  // fall back when configured is missing or only whitespace
  return (trimmed.length > 0 ? trimmed : fallback).replace(/\/+$/, ``)
}

export const join_asset_url = (base_url: string, asset: string): string =>
  `${base_url}/${asset.replace(/^\/+/, ``)}`

// Shared manifest-based URL helpers for parity asset modules (energy, kappa, ...).
// `kind` only affects error messages. `env_base_url` (e.g. a Vite env var) overrides
// the manifest's local default base URL.
export function parity_asset_resolver(
  kind: string,
  manifest: { local_asset_base_url: string; model_assets: object },
  env_base_url: string | undefined,
): {
  asset_url: (asset: string) => string
  model_asset: (model_key: string) => string
  has_model: (model_key: string | undefined) => boolean
} {
  const model_assets = manifest.model_assets as Record<
    string,
    { asset: string } | undefined
  >
  const base_url = resolve_asset_base_url(env_base_url, manifest.local_asset_base_url)
  return {
    asset_url: (asset: string): string => join_asset_url(base_url, asset),
    model_asset: (model_key: string): string => {
      const asset = model_assets[model_key]?.asset
      if (!asset) throw new Error(`No ${kind} parity model asset for ${model_key}`)
      return asset
    },
    has_model: (model_key: string | undefined): boolean =>
      Boolean(model_key && model_assets[model_key]),
  }
}

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

// Keys of TModel whose values are arrays, e.g. the per-row prediction columns. Excludes
// scalar fields like model_key/model_label so only row-aligned fields can be validated.
type ArrayKeys<TModel> = {
  [K in keyof TModel]: TModel[K] extends readonly unknown[] ? K : never
}[keyof TModel] &
  string

// Load a per-model parity asset and validate its model key and the row count of its
// primary prediction column. `kind` (e.g. `energy`/`kappa`) only affects error messages.
export async function load_parity_model<TModel extends ParityModel>(
  kind: string,
  url: string,
  model_key: string,
  pred_field: ArrayKeys<TModel>,
  row_count: number,
): Promise<TModel> {
  const { model } = await load_json_asset<ModelAsset<TModel>>(url)
  if (!model) throw new Error(`No ${kind} parity model ${model_key} in its asset`)
  if (model.model_key !== model_key) {
    throw new Error(
      `Invalid ${kind} parity model: expected ${model_key}, got ${model.model_key}`,
    )
  }
  assert_array_length(
    `${kind} parity ${model_key}.${pred_field}`,
    model[pred_field],
    row_count,
  )
  return model
}
