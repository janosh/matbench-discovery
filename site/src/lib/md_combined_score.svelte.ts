import { MD_METRICS } from '$lib/labels'
import type { Label, ModelData } from '$lib/types'

// Combined MD score (CMDS): like CPS, computed on the fly from the stored component
// errors instead of being persisted in model YAMLs, so it can never go stale when the
// formula changes and users can reweight the components live (e.g. via RadarChart).
// RDF is excluded: 0.9+ cross-model correlation with vDOS/ADF would double-count
// structural accuracy (see matbench_discovery/metrics/md-metrics-design.md).
export const DEFAULT_CMDS_CONFIG = {
  adf_error: { ...MD_METRICS.md_adf_error, weight: 1 / 3 },
  vdos_error: { ...MD_METRICS.md_vdos_error, weight: 1 / 3 },
  pressure_error: { ...MD_METRICS.md_pressure_error, weight: 1 / 3 },
} as const

export type CmdsConfig = Record<
  keyof typeof DEFAULT_CMDS_CONFIG,
  Label & { weight: number }
>
export const CMDS_CONFIG: CmdsConfig = $state(structuredClone(DEFAULT_CMDS_CONFIG))

// Each component is a [0,100] % error (lower=better); its subscore is 1 - error/100,
// clamped to [0,1]. CMDS is the weighted mean of the subscores.
export function calculate_cmds(
  errors: Partial<Record<keyof CmdsConfig, number | undefined>>,
  cmds_config: CmdsConfig,
): number | null {
  const keys = Object.keys(cmds_config) as (keyof CmdsConfig)[]
  const total_weight = keys.reduce((sum, key) => sum + cmds_config[key].weight, 0)
  if (total_weight === 0) return 0

  let weighted_sum = 0
  for (const key of keys) {
    const { weight } = cmds_config[key]
    if (weight === 0) continue
    const error = errors[key]
    // any missing/NaN component with non-zero weight invalidates the score
    if (typeof error !== `number` || isNaN(error)) return null
    weighted_sum += Math.max(0, 1 - error / 100) * weight
  }
  return weighted_sum / total_weight
}

// Write freshly computed CMDS into each model's metrics.md.combined_score so all
// downstream consumers (tables, scatter plots, URL params) keep reading the same path.
export function update_models_cmds(models: ModelData[], cmds_config: CmdsConfig) {
  for (const model of models) {
    const md = model.metrics?.md
    if (!md || typeof md !== `object`) continue
    const { adf_error, vdos_error, pressure_error } = md
    const cmds = calculate_cmds({ adf_error, vdos_error, pressure_error }, cmds_config)
    const md_scored = md as { combined_score?: number }
    // delete rather than write NaN when components are missing: a NaN field still
    // counts as "present" in scatter model counts and exports literally in CSVs
    if (cmds === null) delete md_scored.combined_score
    else md_scored.combined_score = cmds
  }
}
