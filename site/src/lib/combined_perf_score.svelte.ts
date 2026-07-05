import { ALL_METRICS, RMSD_BASELINE } from '$lib/labels'
import type { Label } from '$lib/types'

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

// F1 score is between 0-1 where higher is better (no normalization needed)
function normalize_f1(value: number | undefined): number {
  if (value === undefined || isNaN(value)) return 0
  return value // Already in [0,1] range
}

// RMSD is lower=better, with current models in the range of ~0.01-0.25 (unitless)
// We invert this so that better performance = higher score
function normalize_rmsd(value: number | undefined): number {
  if (value === undefined || isNaN(value)) return 0

  // Fixed reference points for StructureMatcher normalized RMSD:
  // RMS displacement divided by (V/N)^(1/3) (cell volume per atom), so unitless
  const excellent = 0 // Perfect performance (atoms in exact correct positions)

  // Linear interpolation between fixed points with clamping
  // Inverse mapping since lower RMSD is better
  if (value <= excellent) return 1
  if (value >= RMSD_BASELINE) return 0
  return (RMSD_BASELINE - value) / (RMSD_BASELINE - excellent)
}

// κ_SRME is symmetric relative mean error, with range [0,2] by definition
// Lower values are better (0 is perfect)
function normalize_kappa_srme(value: number | undefined): number {
  if (value === undefined || isNaN(value)) return 0

  // Simple linear normalization from [0,2] to [1,0]
  // No clamping needed as SRME is bounded by definition
  return Math.max(0, 1 - value / 2)
}

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
    (F1.weight > 0 && (f1 === undefined || isNaN(f1))) ||
    (RMSD.weight > 0 && (rmsd === undefined || isNaN(rmsd))) ||
    (κ_SRME.weight > 0 && (kappa === undefined || isNaN(kappa)))
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
