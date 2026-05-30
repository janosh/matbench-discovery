import type { AnyStructure } from 'matterviz/structure'
import { parse_any_structure } from 'matterviz/structure/parse'
import {
  assert_array_length,
  clear_asset_cache,
  join_asset_url,
  load_json_asset,
  resolve_asset_base_url,
} from './asset-loader'
import type { ModelAsset, ParityBase, ParityModel, ParityPoint } from './asset-loader'
import { energy_parity_manifest } from './energy-parity-manifest'
import { is_finite_num } from './metrics'

export { load_json_asset }

export type EnergyKind = `e-form` | `each`

export interface EnergyParityBase extends ParityBase {
  n_sites: (number | null)[]
  e_form_true: (number | null)[]
  each_true: (number | null)[]
  structure_shard_size: number
}

export interface EnergyParityModel extends ParityModel {
  e_form_pred: (number | null)[]
}

export interface EnergyParityPoint extends ParityPoint {
  row_idx: number
  n_sites: number
  e_form_true: number
  e_form_pred: number
  each_true: number
  each_pred: number
  x: number
  y: number
  error: number
}

export interface EnergyParitySeries {
  x: Float32Array
  y: Float32Array
  point_ids: Uint32Array
  size_values: Float32Array
}

export type StructurePopupSide = `left` | `right` | `overlay`

export interface StructurePopupPlacement {
  side: StructurePopupSide
  place_right: boolean
  left: number
  top: number
}

export interface StructurePopupPlacementOptions {
  viewport_width: number
  plot_left: number
  plot_width: number
  plot_height: number
  plot_pad_left?: number
  plot_pad_right?: number
  popup_width?: number
  gap?: number
  viewport_margin?: number
}

interface EnergyParityStructureBundle {
  shards: Record<string, Record<string, AnyStructure | string> | undefined>
}

const model_assets = energy_parity_manifest.model_assets as Record<
  string,
  { asset: string } | undefined
>

export const energy_parity_asset_base_url = resolve_asset_base_url(
  import.meta.env.VITE_ENERGY_PARITY_ASSET_BASE_URL as string | undefined,
  energy_parity_manifest.local_asset_base_url,
)

// clears the shared asset cache for all asset types (energy-parity, kappa-parity, ...),
// not just energy-parity; see clear_asset_cache in asset-loader.ts
export const clear_energy_parity_asset_cache = clear_asset_cache

export const energy_parity_asset_url = (asset: string): string =>
  join_asset_url(energy_parity_asset_base_url, asset)

export function structure_popup_placement({
  viewport_width,
  plot_left,
  plot_width,
  plot_height,
  plot_pad_left = 64,
  plot_pad_right = 24,
  popup_width = 500,
  gap = 64,
  viewport_margin = 16,
}: StructurePopupPlacementOptions): StructurePopupPlacement {
  const left_anchor = plot_pad_left
  const right_anchor = Math.max(left_anchor, plot_width - plot_pad_right)
  const required_space = popup_width + gap
  const left_space = plot_left + left_anchor - viewport_margin
  const right_space = viewport_width - (plot_left + right_anchor) - viewport_margin

  const side =
    left_space >= required_space
      ? `left`
      : right_space >= required_space
        ? `right`
        : `overlay`

  const left = side === `left` ? left_anchor : right_anchor
  return { side, place_right: side === `right`, left, top: plot_height / 2 }
}

function assert_energy_parity_base(base: EnergyParityBase): EnergyParityBase {
  const { row_count, structure_shard_size } = energy_parity_manifest
  for (const key of [
    `material_ids`,
    `formulas`,
    `n_sites`,
    `e_form_true`,
    `each_true`,
  ] as const) {
    assert_array_length(`energy parity ${key}`, base[key], row_count)
  }
  if (base.structure_shard_size !== structure_shard_size) {
    throw new Error(
      `Invalid energy parity structure shard size: expected ${structure_shard_size}, ` +
        `got ${base.structure_shard_size}`,
    )
  }
  return base
}

function assert_energy_parity_model(
  model: EnergyParityModel,
  model_key: string,
): EnergyParityModel {
  if (model.model_key !== model_key) {
    throw new Error(
      `Invalid energy parity model: expected ${model_key}, got ${model.model_key}`,
    )
  }
  assert_array_length(
    `energy parity ${model_key}.e_form_pred`,
    model.e_form_pred,
    energy_parity_manifest.row_count,
  )
  return model
}

export function model_asset(model_key: string): string {
  const asset = model_assets[model_key]?.asset
  if (!asset) throw new Error(`No energy parity model asset for ${model_key}`)
  return asset
}

export const has_energy_parity_model = (
  model: EnergyParityModel | undefined,
  model_key: string | undefined,
): model is EnergyParityModel => Boolean(model_key && model?.model_key === model_key)

export const load_energy_parity_base = (): Promise<EnergyParityBase> =>
  load_json_asset<EnergyParityBase>(
    energy_parity_asset_url(energy_parity_manifest.base.asset),
  ).then(assert_energy_parity_base)

export async function load_energy_parity_model(
  model_key: string,
): Promise<EnergyParityModel> {
  const { model } = await load_json_asset<ModelAsset<EnergyParityModel>>(
    energy_parity_asset_url(model_asset(model_key)),
  )
  if (!model) throw new Error(`No energy parity model ${model_key} in its asset`)
  return assert_energy_parity_model(model, model_key)
}

export function get_energy_parity_point(
  base: EnergyParityBase,
  model: EnergyParityModel,
  row_idx: number,
  energy_kind: EnergyKind,
): EnergyParityPoint | null {
  const e_form_true = base.e_form_true[row_idx]
  const each_true = base.each_true[row_idx]
  const e_form_pred = model.e_form_pred[row_idx]
  const n_sites = base.n_sites[row_idx]
  if (
    !is_finite_num(e_form_true) ||
    !is_finite_num(each_true) ||
    !is_finite_num(e_form_pred) ||
    !is_finite_num(n_sites)
  ) {
    return null
  }

  const each_pred = each_true + e_form_pred - e_form_true
  const is_each = energy_kind === `each`
  const x = is_each ? each_true : e_form_true
  const y = is_each ? each_pred : e_form_pred
  if (!is_finite_num(x) || !is_finite_num(y)) return null

  return {
    row_idx,
    material_id: base.material_ids[row_idx] ?? `unknown-${row_idx}`,
    formula: base.formulas[row_idx] ?? ``,
    n_sites,
    e_form_true,
    e_form_pred,
    each_true,
    each_pred,
    x,
    y,
    error: y - x,
  }
}

export function build_energy_parity_series(
  base: EnergyParityBase,
  model: EnergyParityModel,
  energy_kind: EnergyKind,
): EnergyParitySeries {
  // get_energy_parity_point returns null for any out-of-range/non-finite row,
  // so iterating material_ids covers all valid rows regardless of array lengths
  const points = base.material_ids
    .map((_material_id, row_idx) =>
      get_energy_parity_point(base, model, row_idx, energy_kind),
    )
    .filter((point): point is EnergyParityPoint => point !== null)

  return {
    x: Float32Array.from(points, (point) => point.x),
    y: Float32Array.from(points, (point) => point.y),
    point_ids: Uint32Array.from(points, (point) => point.row_idx),
    size_values: Float32Array.from(points, (point) => point.n_sites),
  }
}

export function structure_shard_idx(
  base: Pick<EnergyParityBase, `structure_shard_size`>,
  row_idx: number,
): number {
  if (!Number.isInteger(row_idx) || row_idx < 0)
    throw new Error(`Invalid row index ${row_idx}`)
  if (!Number.isInteger(base.structure_shard_size) || base.structure_shard_size <= 0) {
    throw new Error(`Invalid structure shard size ${base.structure_shard_size}`)
  }
  return Math.floor(row_idx / base.structure_shard_size)
}

export function structure_bundle_for_shard(shard_idx: number) {
  const bundle = energy_parity_manifest.structure_bundles.find(
    ({ start_shard, end_shard }) => start_shard <= shard_idx && shard_idx <= end_shard,
  )
  if (!bundle) throw new Error(`No structure bundle for shard ${shard_idx}`)
  return bundle
}

export async function load_wbm_structure(
  base: Pick<EnergyParityBase, `structure_shard_size`>,
  row_idx: number,
  material_id: string,
): Promise<AnyStructure> {
  const shard_idx = structure_shard_idx(base, row_idx)
  const bundle = await load_json_asset<EnergyParityStructureBundle>(
    energy_parity_asset_url(structure_bundle_for_shard(shard_idx).asset),
  )
  const payload = bundle.shards?.[String(shard_idx)]?.[material_id]
  if (!payload) throw new Error(`No structure found for ${material_id}`)
  if (typeof payload !== `string`) return payload

  const structure = parse_any_structure(payload, `${material_id}.extxyz`)
  if (!structure) throw new Error(`Could not parse structure for ${material_id}`)
  return structure
}
