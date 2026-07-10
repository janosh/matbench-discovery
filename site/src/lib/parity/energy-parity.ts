import { fsum, variance } from 'd3-array'
import type { AnyStructure } from 'matterviz/structure'
import { parse_any_structure } from 'matterviz/structure/parse'
import {
  assert_array_length,
  load_json_asset,
  load_parity_model,
  parity_asset_resolver,
} from '../asset-loader'
import type { ParityBase, ParityModel, ParityPoint } from '../asset-loader'
import { energy_parity_manifest } from './energy-parity-manifest'
import { is_finite_num } from '../metrics'

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

export interface StructurePopupPlacement {
  side: `left` | `right`
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

export const { asset_url: energy_parity_asset_url, model_asset } = parity_asset_resolver(
  `energy`,
  energy_parity_manifest,
  import.meta.env.VITE_ENERGY_PARITY_ASSET_BASE_URL,
)

export function structure_popup_placement({
  viewport_width,
  plot_left,
  plot_width,
  plot_height,
  plot_pad_left = 64,
  plot_pad_right = 24,
  popup_width = 500,
  gap = 16,
  viewport_margin = 16,
}: StructurePopupPlacementOptions): StructurePopupPlacement {
  // anchors sit at the plot's data-area edges (the axis padding has no points)
  const left_anchor = plot_pad_left
  const right_anchor = Math.max(left_anchor, plot_width - plot_pad_right)
  // pick the roomier side so the popup overlaps as little of the plot as possible;
  // ties go left (over the y axis). Never centers over the data anymore.
  const left_space = plot_left + left_anchor
  const right_space = viewport_width - (plot_left + right_anchor)
  const place_right = right_space > left_space

  let left: number
  if (place_right) {
    // popup hangs right of the anchor (left edge at plot_left + anchor + gap); slide it
    // out to the data edge, but clamp left so its right edge stays in the viewport
    const max_anchor = viewport_width - viewport_margin - popup_width - gap - plot_left
    left = Math.min(right_anchor, max_anchor)
  } else {
    // popup hangs left of the anchor (right edge at plot_left + anchor - gap); slide it
    // out to the data edge, but clamp right so its left edge stays in the viewport
    const min_anchor = viewport_margin - plot_left + gap + popup_width
    left = Math.max(left_anchor, min_anchor)
  }

  return { side: place_right ? `right` : `left`, left, top: plot_height / 2 }
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

export const load_energy_parity_base = (): Promise<EnergyParityBase> =>
  load_json_asset<EnergyParityBase>(
    energy_parity_asset_url(energy_parity_manifest.base.asset),
  ).then(assert_energy_parity_base)

export const load_energy_parity_model = (model_key: string): Promise<EnergyParityModel> =>
  load_parity_model<EnergyParityModel>(
    `energy`,
    energy_parity_asset_url(model_asset(model_key)),
    model_key,
    `e_form_pred`,
    energy_parity_manifest.row_count,
  )

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
  const row_count = base.material_ids.length
  const x = new Float32Array(row_count)
  const y = new Float32Array(row_count)
  const point_ids = new Uint32Array(row_count)
  const size_values = new Float32Array(row_count)
  let valid_count = 0

  for (let row_idx = 0; row_idx < row_count; row_idx++) {
    const point = get_energy_parity_point(base, model, row_idx, energy_kind)
    if (!point) continue
    x[valid_count] = point.x
    y[valid_count] = point.y
    point_ids[valid_count] = row_idx
    size_values[valid_count] = point.n_sites
    valid_count++
  }
  return {
    x: x.subarray(0, valid_count),
    y: y.subarray(0, valid_count),
    point_ids: point_ids.subarray(0, valid_count),
    size_values: size_values.subarray(0, valid_count),
  }
}

export interface EnergyParityStats {
  mae: number // mean absolute error of ML vs DFT energy
  r2: number // coefficient of determination, R^2 = 1 - SS_res / SS_tot
  n_points: number
}

// MAE and R^2 of the plotted DFT (x) vs ML (y) energies, for the in-plot annotation.
// Computed from the displayed points so it matches whatever energy_kind is shown.
export function energy_parity_stats(series: EnergyParitySeries): EnergyParityStats {
  const { x, y } = series
  const n_points = x.length
  if (n_points === 0) return { mae: NaN, r2: NaN, n_points: 0 }

  // fsum = compensated (full-precision) summation, so the aggregate is exact over the
  // ~256k terms; variance (Welford) gives a stable SS_tot without a separate mean pass
  const mae = fsum(x, (x_val, idx) => Math.abs(y[idx] - x_val)) / n_points
  const ss_res = fsum(x, (x_val, idx) => (y[idx] - x_val) ** 2)
  const ss_tot = (variance(x) ?? 0) * (n_points - 1) // d3 variance is sample (n-1)
  return { mae, r2: ss_tot > 0 ? 1 - ss_res / ss_tot : NaN, n_points }
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
