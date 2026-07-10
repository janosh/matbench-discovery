import { ALL_METRICS, MD_METRICS, RMSD_BASELINE } from '$lib/labels'
import type { DiatomicsMetrics, MdMetrics } from '$lib/schema/model'
import type { Label, ModelData } from '$lib/types'

export const DEFAULT_CPS_CONFIG = {
  F1: { ...ALL_METRICS.F1, weight: 0.5 },
  κ_SRME: { ...ALL_METRICS.κ_SRME, weight: 0.4 },
  RMSD: { ...ALL_METRICS.RMSD, weight: 0.1 },
} as const

export type CpsConfig = Record<
  keyof typeof DEFAULT_CPS_CONFIG,
  Label & { weight: number }
>
// deep clone: a shallow spread would share the nested metric objects with
// DEFAULT_CPS_CONFIG, coupling weight edits to the defaults they're reset from
export const CPS_CONFIG: CpsConfig = $state(structuredClone(DEFAULT_CPS_CONFIG))

const is_valid_score = (value: number | undefined): value is number =>
  value !== undefined && !isNaN(value)

// F1 score is between 0-1 where higher is better (no normalization needed)
const normalize_f1 = (value: number | undefined): number =>
  is_valid_score(value) ? value : 0

// RMSD is lower=better, with current models in the range of ~0.01-0.25 (unitless)
// We invert this so that better performance = higher score
const normalize_rmsd = (value: number | undefined): number =>
  is_valid_score(value) ? Math.max(0, Math.min(1, 1 - value / RMSD_BASELINE)) : 0

// κ_SRME is symmetric relative mean error, with range [0,2] by definition
// Lower values are better (0 is perfect)
const normalize_kappa_srme = (value: number | undefined): number =>
  is_valid_score(value) ? Math.max(0, 1 - value / 2) : 0

// Calculate a combined score using normalized metrics weighted by importance factors.
// Fixed normalization reference points keep scores stable as new models are added.
export function calculate_cps(
  f1: number | undefined,
  rmsd: number | undefined,
  kappa: number | undefined,
  cps_config: CpsConfig, // Weights for each metric
): number | null {
  const { F1, RMSD, κ_SRME } = cps_config

  // Any metric with non-zero weight must be present, else the score is undefined
  if (
    (F1.weight > 0 && !is_valid_score(f1)) ||
    (RMSD.weight > 0 && !is_valid_score(rmsd)) ||
    (κ_SRME.weight > 0 && !is_valid_score(kappa))
  )
    return null

  // all-zero weights leave the score undefined rather than tying every model at 0
  // (matches calculate_cmds; unreachable via UI, which rejects all-zero weights)
  const total_weight = F1.weight + RMSD.weight + κ_SRME.weight
  if (total_weight === 0) return null

  // zero-weight metrics contribute 0 regardless of value (normalizers map undefined to 0)
  const weighted_sum =
    normalize_f1(f1) * F1.weight +
    normalize_rmsd(rmsd) * RMSD.weight +
    normalize_kappa_srme(kappa) * κ_SRME.weight

  return weighted_sum / total_weight
}

// Combined MD score (CMDS): like CPS/CDS, computed on the fly from the stored
// components instead of being persisted in model YAMLs, so it can never go stale when
// the formula changes and users can reweight the components live (e.g. via RadarChart).
// RDF is excluded: 0.9+ cross-model correlation with vDOS/ADF would double-count
// structural accuracy (see matbench_discovery/metrics/md-metrics-design.md).
//
// Speed = summed NVT rollout wall time, log-clamped from floor (subscore 1) to
// baseline (subscore 0) like the CDS speed pillar. Bounds are FROZEN just under the
// fastest model (~9.4e3 s) and near the 90th percentile (~2.6e5 s) of the mid-2026
// field (33 models, one NVIDIA H200 per system) so scores stay stable as models are
// added.
export const CMDS_SPEED = { baseline: 300_000, floor: 9000, log: true } as const

// Key order = RadarChart corner order (cyclic; diagonals vDOS↔speed, ADF↔pressure).
// Default weights must be knob-expressible (see DEFAULT_CDS_CONFIG): opposite-corner
// products must match, which 30/20/20/30 satisfies (0.3·0.2 == 0.2·0.3). vDOS and
// pressure get the emphasis: vDOS is the headline dynamical observable and pressure
// the least correlated with the other components (see md-metrics-design.md).
export const DEFAULT_CMDS_CONFIG = {
  vdos_error: { ...MD_METRICS.md_vdos_error, weight: 0.3 },
  adf_error: { ...MD_METRICS.md_adf_error, weight: 0.2 },
  run_time_sec: {
    ...MD_METRICS.md_run_time_sec,
    label: `Speed`,
    description: `Wall time to roll out all 17 DynaMat v1.0 NVT trajectories, scored on a log scale from 9,000 s (subscore 1) to 300,000 s (subscore 0). All timings to date were measured on one NVIDIA H200 per system; models without recorded timings get no CMDS unless this weight is zeroed`,
    weight: 0.2,
  },
  pressure_error: { ...MD_METRICS.md_pressure_error, weight: 0.3 },
} as const

export type CmdsConfig = Record<
  keyof typeof DEFAULT_CMDS_CONFIG,
  Label & { weight: number }
>
export const CMDS_CONFIG: CmdsConfig = $state(structuredClone(DEFAULT_CMDS_CONFIG))

// Clamped [0,1] subscore for a lower=better value on a floor->baseline scale,
// optionally in log10 space. Shared by the CMDS speed component and all CDS pillars.
type ScoreScale = { baseline: number; floor: number; log: boolean }
const subscore = ({ baseline, floor, log }: ScoreScale, value: number): number => {
  const [val, lo, hi] = log
    ? [Math.log10(value), Math.log10(floor), Math.log10(baseline)]
    : [value, floor, baseline]
  return Math.min(1, Math.max(0, 1 - (val - lo) / (hi - lo)))
}

// Error components are [0,100] % errors (lower=better) with subscore 1 - error/100;
// run_time_sec is scored on the CMDS_SPEED log scale. CMDS is the weighted mean of
// the subscores.
export function calculate_cmds(
  values: Partial<Record<keyof CmdsConfig, number>>,
  cmds_config: CmdsConfig,
): number | null {
  const keys = Object.keys(cmds_config) as (keyof CmdsConfig)[]
  const total_weight = keys.reduce((sum, key) => sum + cmds_config[key].weight, 0)
  // all-zero weights leave the score undefined, not 0 (which would rank as worst)
  if (total_weight === 0) return null

  let weighted_sum = 0
  for (const key of keys) {
    const { weight } = cmds_config[key]
    if (weight === 0) continue
    const value = values[key]
    // any missing/invalid component with non-zero weight invalidates the score
    if (typeof value !== `number` || !isFinite(value)) return null
    if (key === `run_time_sec`) {
      if (value <= 0) return null // log scale needs a positive wall time
      weighted_sum += subscore(CMDS_SPEED, value) * weight
    } else weighted_sum += Math.max(0, 1 - value / 100) * weight
  }
  return weighted_sum / total_weight
}

// Recompute a task's on-the-fly combined_score into every model's
// metrics[task].combined_score, the shared path all consumers (tables, scatter plots,
// URL params) read. Delete rather than write NaN when components are missing: a NaN
// field still counts as "present" in scatter model counts and exports literally in CSVs.
function write_combined_scores(
  models: ModelData[],
  task: `md` | `diatomics`,
  score: (metrics: DiatomicsMetrics | MdMetrics) => number | null,
) {
  for (const model of models) {
    const metrics = model.metrics?.[task]
    if (!metrics || typeof metrics !== `object`) continue
    const value = score(metrics)
    if (value === null) Reflect.deleteProperty(metrics, `combined_score`)
    else Object.assign(metrics, { combined_score: value })
  }
}

export const update_models_cmds = (models: ModelData[], config: CmdsConfig) =>
  write_combined_scores(models, `md`, (metrics) => calculate_cmds(metrics, config))

// Combined Diatomics Score (CDS): like CPS/CMDS, computed on the fly from the stored
// per-metric errors instead of being persisted in model YAMLs, so it can never go
// stale when the formula changes and users can reweight the pillars live.
//
// Each component subscore = clamp(1 - (value - floor) / (baseline - floor), 0, 1),
// with value/floor/baseline in log10 space for log: true components (run time spans
// two orders of magnitude, 57s-4500s across models, so a linear clamp would compress
// most of the field into the top decile; log treats 2x slower the same at any scale).
// Baselines are frozen reference scales so scores stay stable as models are added.
// floor=1 for metrics with a physical minimum of 1: tortuosity (path length /
// net displacement) and force_flips (a well-formed PEC has exactly one force
// zero-crossing at the equilibrium bond length).
//
// Deliberately excluded: pbe_vib_freq_error (curvature at the minimum, already
// captured by force MAE + well depth, heavy-tailed 43-450 cm⁻¹ across models),
// energy_diff_flips (near-duplicate of force_flips), force_total_variation and
// force_jump (double-count smoothness already scored by energy_jump/force_flips and
// can penalize physically sharp repulsive walls).
// component weights sum to 1 within each pillar. Pillar (= RadarChart corner) order
// is cyclic: accuracy and speed sit on opposite corners, as do geometry and
// physicality - see DEFAULT_CDS_CONFIG for why that pairing matters
export const CDS_COMPONENTS = {
  accuracy: [
    { key: `pbe_energy_mae`, weight: 5 / 9, baseline: 4, floor: 0, log: false }, // eV
    { key: `pbe_force_mae`, weight: 4 / 9, baseline: 5, floor: 0, log: false }, // eV/Å
  ],
  geometry: [
    { key: `pbe_wall_dist_mae`, weight: 1 / 2, baseline: 0.5, floor: 0, log: false }, // Å
    { key: `pbe_bond_length_error`, weight: 1 / 3, baseline: 0.8, floor: 0, log: false }, // Å
    { key: `pbe_well_depth_error`, weight: 1 / 6, baseline: 3.5, floor: 0, log: false }, // eV
  ],
  speed: [
    { key: `run_time_sec`, weight: 1, baseline: 5000, floor: 50, log: true }, // s
  ],
  physicality: [
    { key: `energy_jump`, weight: 2 / 5, baseline: 3, floor: 0, log: false }, // eV
    { key: `force_flips`, weight: 2 / 5, baseline: 3.5, floor: 1, log: false },
    { key: `tortuosity`, weight: 1 / 5, baseline: 2, floor: 1, log: false },
  ],
} as const

export type CdsPillar = keyof typeof CDS_COMPONENTS
export type CdsConfig = Record<CdsPillar, Label & { weight: number }>
type CdsComponent = (typeof CDS_COMPONENTS)[CdsPillar][number]
export type CdsValues = Partial<Record<CdsComponent[`key`], number>> & {
  excluded_formula_reasons?: Record<string, string>
}
const N_SCORED_DIATOMIC_ELEMENTS = 87

// Default weights MUST be a weight vector the RadarChart knob can express, else
// Reset shows a knob position whose canonical reading disagrees with the displayed
// percentages. A 2D knob has 2 degrees of freedom vs 3 for four weights: on a square
// the reachable (Wachspress = bilinear) weights are exactly the products of two
// marginal splits, so opposite-corner products must match:
// w_accuracy·w_speed == w_geometry·w_physicality. 4/9·1/9 == 2/9·2/9 satisfies this
// (marginal splits 2/3 toward accuracy on both axes), approximating the intended
// 40/25/25/10 while keeping Reset, knob position and weights mutually consistent.
export const DEFAULT_CDS_CONFIG: CdsConfig = {
  accuracy: {
    key: `cds_accuracy`,
    label: `Accuracy`,
    description: `PEC agreement with PBE: energy MAE (5/9) and force MAE (4/9) over the scored separation range. Note PBE itself is unreliable for stretched diatomics (dissociation limits, spin states), so models trained on other references may deviate here without being worse potentials`,
    weight: 4 / 9,
  },
  geometry: {
    key: `cds_geometry`,
    label: `Geometry`,
    description: `Curve-shape agreement with PBE: repulsive-wall distance from 1 to 100 eV where supported by the reference (1/2), equilibrium bond length (1/3) and well depth (1/6) errors`,
    weight: 2 / 9,
  },
  speed: {
    key: `cds_speed`,
    label: `Speed`,
    description: `Wall time of the full diatomic sweep, scored on a log scale from 50s (subscore 1) to 5000s (subscore 0). Caveat: hardware varies by submission (recorded in each model YAML); models without recorded timings get no CDS unless this weight is zeroed`,
    weight: 1 / 9,
  },
  physicality: {
    key: `cds_physicality`,
    label: `Physicality`,
    description: `Reference-free smoothness diagnostics: energy jumps at sign flips (2/5), force-direction flips (2/5) and tortuosity (1/5). Penalizes pathological potential energy curves regardless of DFT reference`,
    weight: 2 / 9,
  },
}

// deep clone so weight edits don't mutate the defaults they're reset from
export const CDS_CONFIG: CdsConfig = $state(structuredClone(DEFAULT_CDS_CONFIG))

// CDS is the weighted mean of the pillar subscores, each pillar a fixed-weight mean
// of its component subscores, multiplied by the fraction of the 87 benchmark elements
// the model completed. Any missing/NaN component in a non-zero-weight pillar invalidates
// the score (matches calculate_cps/calculate_cmds).
export function calculate_cds(values: CdsValues, cds_config: CdsConfig): number | null {
  const pillars = Object.keys(cds_config) as CdsPillar[]
  const total_weight = pillars.reduce((sum, key) => sum + cds_config[key].weight, 0)
  // all-zero weights leave the score undefined, not 0 (which would rank as worst)
  if (total_weight === 0) return null

  let weighted_sum = 0
  for (const pillar of pillars) {
    const { weight } = cds_config[pillar]
    if (weight === 0) continue
    let pillar_score = 0
    for (const component of CDS_COMPONENTS[pillar]) {
      const value = values[component.key]
      if (typeof value !== `number` || !isFinite(value) || (component.log && value <= 0))
        return null
      pillar_score += subscore(component, value) * component.weight
    }
    weighted_sum += pillar_score * weight
  }
  const n_excluded = Object.keys(values.excluded_formula_reasons ?? {}).length
  const coverage = Math.max(0, 1 - n_excluded / N_SCORED_DIATOMIC_ELEMENTS)
  return (weighted_sum / total_weight) * coverage
}

export const update_models_cds = (models: ModelData[], config: CdsConfig) =>
  write_combined_scores(models, `diatomics`, (metrics) => calculate_cds(metrics, config))
