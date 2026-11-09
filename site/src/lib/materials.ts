// Per-material views of the κ-103 analysis payload: one row per model for a single
// phononDB material, plus cross-model summaries for the /materials index. Pure
// functions so prerender loaders can slice the 1 MB payload down to what a page needs.
import type { KappaAnalysis } from '#lib/parity/kappa-parity.js'
import { is_finite_num } from '#lib/metrics.js'
import { SPACEGROUP_NUM_TO_SYMBOL } from 'matterviz/symmetry'

// Hermann-Mauguin symbol as HTML: screw-axis subscripts (P6_3mc -> P6<sub>3</sub>mc)
// and rotoinversion overbars (Fm-3m -> Fm3̄m), followed by the spacegroup number
export function spacegroup_html(spg_num: number): string {
  const symbol = SPACEGROUP_NUM_TO_SYMBOL[spg_num]
  if (!symbol) return String(spg_num)
  const html = symbol
    .replaceAll(/_(?<digit>\d)/g, `<sub>$<digit></sub>`)
    .replaceAll(
      /-(?<digit>\d)/g,
      `<span style="text-decoration: overline">$<digit></span>`,
    )
  return `${html} (${spg_num})`
}

export type KappaMaterialModelRow = {
  model_key: string
  label: string
  kappa_ml: number | null
  kappa_ratio: number | null // κ_ML / κ_DFT
  srme: number | null
  srme_censored: boolean
  freq_w1: number | null
  imag_modes: boolean
  broken_sym: boolean
  max_steps: boolean
  below_floor: boolean
  kappa_min: number | null
}

export type KappaMaterialSummary = {
  material_id: string
  formula: string
  spg_num: number
  kappa_dft: number | null
  n_models: number // models with any prediction for this material
  n_censored: number // SRME censored to 2 (imaginary modes, broken symmetry, failed κ)
  median_srme: number | null
  median_kappa_ratio: number | null
  best: { model_key: string; label: string; srme: number } | null
  worst: { model_key: string; label: string; srme: number } | null
}

export type KappaMaterial = KappaMaterialSummary & {
  rows: KappaMaterialModelRow[] // sorted by SRME, best first; models without a prediction last
  // 1-based rank of median_srme among all materials with one (1 = lowest median error)
  median_srme_rank: number | null
  n_ranked: number
  prev_id: string | null
  next_id: string | null
}

const median = (values: number[]): number | null => {
  if (!values.length) return null
  const sorted = values.toSorted((val_a, val_b) => val_a - val_b)
  const mid = Math.floor(sorted.length / 2)
  return sorted.length % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2
}

// null SRME (no prediction) sorts after every scored model
const by_srme = (row_a: KappaMaterialModelRow, row_b: KappaMaterialModelRow): number =>
  (row_a.srme ?? Infinity) - (row_b.srme ?? Infinity) ||
  row_a.label.localeCompare(row_b.label)

function model_rows(analysis: KappaAnalysis, idx: number): KappaMaterialModelRow[] {
  const kappa_dft = analysis.kappa_dft[idx]
  return analysis.models
    .map((model) => {
      const kappa_ml = model.kappa_ml[idx]
      return {
        model_key: model.model_key,
        label: model.label,
        kappa_ml,
        kappa_ratio:
          is_finite_num(kappa_ml) && is_finite_num(kappa_dft) && kappa_dft > 0
            ? kappa_ml / kappa_dft
            : null,
        srme: model.srme[idx],
        srme_censored: model.srme_censored[idx] === true,
        freq_w1: model.freq_w1[idx],
        imag_modes: model.imag_modes[idx] === true,
        broken_sym: model.broken_sym[idx] === true,
        max_steps: model.max_steps[idx] === true,
        below_floor: model.below_floor[idx] === true,
        kappa_min: model.kappa_min[idx],
      }
    })
    .toSorted(by_srme)
}

function summarize(
  analysis: KappaAnalysis,
  idx: number,
  rows: KappaMaterialModelRow[],
): KappaMaterialSummary {
  const scored = rows.filter((row) => is_finite_num(row.srme))
  const extreme = (row: KappaMaterialModelRow | undefined) =>
    row && is_finite_num(row.srme)
      ? { model_key: row.model_key, label: row.label, srme: row.srme }
      : null
  return {
    material_id: analysis.material_ids[idx],
    formula: analysis.formulas[idx],
    spg_num: analysis.spg_nums[idx],
    kappa_dft: analysis.kappa_dft[idx],
    n_models: rows.filter((row) => row.srme !== null || row.kappa_ml !== null).length,
    n_censored: rows.filter((row) => row.srme_censored).length,
    median_srme: median(scored.flatMap((row) => row.srme ?? [])),
    median_kappa_ratio: median(rows.flatMap((row) => row.kappa_ratio ?? [])),
    best: extreme(scored[0]),
    worst: extreme(scored.at(-1)),
  }
}

function material_idx(analysis: KappaAnalysis, material_id: string): number {
  const idx = analysis.material_ids.indexOf(material_id)
  if (idx === -1) throw new Error(`Material ${material_id} not in the κ-103 analysis`)
  return idx
}

export function kappa_material(
  analysis: KappaAnalysis,
  material_id: string,
): KappaMaterial {
  const idx = material_idx(analysis, material_id)
  const rows = model_rows(analysis, idx)
  const { material_ids } = analysis
  const summary = summarize(analysis, idx, rows)
  const medians = kappa_material_summaries(analysis).flatMap(
    ({ median_srme }) => median_srme ?? [],
  )
  const { median_srme } = summary
  return {
    ...summary,
    rows,
    median_srme_rank:
      median_srme === null
        ? null
        : medians.filter((other) => other < median_srme).length + 1,
    n_ranked: medians.length,
    prev_id: material_ids[idx - 1] ?? null,
    next_id: material_ids[idx + 1] ?? null,
  }
}

export const kappa_material_summaries = (
  analysis: KappaAnalysis,
): KappaMaterialSummary[] =>
  analysis.material_ids.map((_material_id, idx) =>
    summarize(analysis, idx, model_rows(analysis, idx)),
  )
