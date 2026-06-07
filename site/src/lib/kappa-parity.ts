import type { AnyStructure } from 'matterviz/structure'
import { parse_any_structure } from 'matterviz/structure/parse'
import {
  assert_array_length,
  load_parity_model,
  clear_asset_cache,
  join_asset_url,
  load_json_asset,
  resolve_asset_base_url,
} from './asset-loader'
import type { ParityBase, ParityModel, ParityPoint } from './asset-loader'
import { kappa_parity_manifest } from './kappa-parity-manifest'
import { is_finite_num } from './metrics'

export { clear_asset_cache as clear_kappa_parity_asset_cache }

// raw phonon DOS as stored in assets (histogram of mesh frequencies in THz)
export interface RawDos {
  frequencies: number[]
  densities: number[]
}

// shape matterviz `Dos` expects for a phonon density of states
export interface PhononDos extends RawDos {
  type: `phonon`
}

export interface KappaParityBase extends ParityBase {
  kappa_dft: (number | null)[]
  // optional metadata for point styling (size by n_sites, color by crystal system);
  // older assets may lack these, so the client falls back to default size/color
  n_sites?: (number | null)[]
  spacegroups?: (number | null)[]
  structures: Record<string, AnyStructure | string | undefined>
  dft_dos: Record<string, RawDos | undefined>
}

export interface KappaParityModel extends ParityModel {
  kappa_ml: (number | null)[]
  ml_dos: Record<string, RawDos | undefined>
}

export interface KappaParityPoint extends ParityPoint {
  kappa_dft: number
  kappa_ml: number
  sre: number // symmetric relative error of scalar conductivity
  n_sites: number | null
  spacegroup: number | null
}

export interface KappaParitySeries {
  x: number[]
  y: number[]
  points: KappaParityPoint[]
}

const model_assets = kappa_parity_manifest.model_assets as Record<
  string,
  { asset: string } | undefined
>

const kappa_parity_asset_base_url = resolve_asset_base_url(
  import.meta.env.VITE_KAPPA_PARITY_ASSET_BASE_URL as string | undefined,
  kappa_parity_manifest.local_asset_base_url,
)

export const kappa_parity_asset_url = (asset: string): string =>
  join_asset_url(kappa_parity_asset_base_url, asset)

export const kappa_model_asset = (model_key: string): string => {
  const asset = model_assets[model_key]?.asset
  if (!asset) throw new Error(`No kappa parity model asset for ${model_key}`)
  return asset
}

export const has_kappa_parity_model = (model_key: string | undefined): boolean =>
  Boolean(model_key && model_assets[model_key])

export async function load_kappa_parity_base(): Promise<KappaParityBase> {
  const base = await load_json_asset<KappaParityBase>(
    kappa_parity_asset_url(kappa_parity_manifest.base.asset),
  )
  const { row_count } = kappa_parity_manifest
  for (const key of [`material_ids`, `formulas`, `kappa_dft`] as const) {
    assert_array_length(`kappa parity ${key}`, base[key], row_count)
  }
  return base
}

export const load_kappa_parity_model = (model_key: string): Promise<KappaParityModel> =>
  load_parity_model<KappaParityModel>(
    `kappa`,
    kappa_parity_asset_url(kappa_model_asset(model_key)),
    model_key,
    `kappa_ml`,
    kappa_parity_manifest.row_count,
  )

export function get_kappa_parity_point(
  base: KappaParityBase,
  model: KappaParityModel,
  row_idx: number,
): KappaParityPoint | null {
  const x = base.kappa_dft[row_idx]
  const y = model.kappa_ml[row_idx]
  // require positive conductivities: physically meaningful and log-scale safe
  if (!is_finite_num(x) || !is_finite_num(y) || x <= 0 || y <= 0) return null

  return {
    material_id: base.material_ids[row_idx] ?? `unknown-${row_idx}`,
    formula: base.formulas[row_idx] ?? ``,
    kappa_dft: x,
    kappa_ml: y,
    sre: Math.abs((2 * (y - x)) / (x + y)),
    n_sites: base.n_sites?.[row_idx] ?? null,
    spacegroup: base.spacegroups?.[row_idx] ?? null,
  }
}

export function build_kappa_parity_series(
  base: KappaParityBase,
  model: KappaParityModel,
): KappaParitySeries {
  const points = base.material_ids
    .map((_material_id, row_idx) => get_kappa_parity_point(base, model, row_idx))
    .filter((pt): pt is KappaParityPoint => pt !== null)
  return {
    x: points.map((pt) => pt.kappa_dft),
    y: points.map((pt) => pt.kappa_ml),
    points,
  }
}

function as_phonon_dos(dos: RawDos | undefined): PhononDos | null {
  if (!dos?.frequencies?.length || !dos.densities?.length) return null
  return { type: `phonon`, frequencies: dos.frequencies, densities: dos.densities }
}

export const dft_phonon_dos = (
  base: KappaParityBase,
  material_id: string,
): PhononDos | null => as_phonon_dos(base.dft_dos[material_id])

export const ml_phonon_dos = (
  model: KappaParityModel,
  material_id: string,
): PhononDos | null => as_phonon_dos(model.ml_dos[material_id])

export function kappa_structure(
  base: KappaParityBase,
  material_id: string,
): AnyStructure | null {
  const payload = base.structures[material_id]
  if (!payload) return null
  if (typeof payload !== `string`) return payload

  // runs inside a $derived in the component, so degrade gracefully on parse
  // failure (returning null hides the structure) instead of crashing the plot
  try {
    return parse_any_structure(payload, `${material_id}.extxyz`)
  } catch (error) {
    console.warn(`Failed to parse structure for ${material_id}:`, error)
    return null
  }
}
