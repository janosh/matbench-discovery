import { ALL_METRICS, RMSD_BASELINE } from './labels'
import type { Label } from './types'

export const DEFAULT_CPS_CONFIG = {
  F1: { ...ALL_METRICS.F1, weight: 0.5 },
  κ_SRME: { ...ALL_METRICS.κ_SRME, weight: 0.4 },
  RMSD: { ...ALL_METRICS.RMSD, weight: 0.1 },
} as const

export type CpsConfig = Record<
  keyof typeof DEFAULT_CPS_CONFIG,
  Label & { weight: number }
>
// Make CPS_CONFIG reactive (using Svelte 5 runes)
export const CPS_CONFIG: CpsConfig = $state({ ...DEFAULT_CPS_CONFIG })

// F1 score is between 0-1 where higher is better (no normalization needed)
function normalize_f1(value: number | undefined): number {
  if (value === undefined || isNaN(value)) return 0
  return value // Already in [0,1] range
}

// RMSD is lower=better, with current models in the range of ~0.01-0.25 Å
// We invert this so that better performance = higher score
function normalize_rmsd(value: number | undefined): number {
  if (value === undefined || isNaN(value)) return 0

  // Fixed reference points for RMSD (in Å)
  const excellent = 0 // Perfect performance (atoms in exact correct positions)

  // Linear interpolation between fixed points with clamping
  // Inverse mapping since lower RMSD is better
  if (value <= excellent) return 1.0
  if (value >= RMSD_BASELINE) return 0.0
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
// This uses fixed normalization reference points to ensure score stability when new models are added.

// Normalization reference points:
// - F1 score for discovery already in [0,1] range, higher is better
// - RMSD Root mean square displacement in range 0 (perfect) to RMSD_BASELINE, lower is better
// - κ_SRME symmetric relative mean error for lattice thermal conductivity,
//    range [0,2] linearly mapped to [1,0], lower is better
export function calculate_cps(
  f1: number | undefined,
  rmsd: number | undefined,
  kappa: number | undefined,
  cps_config: CpsConfig, // weights for each metric
): number | null {
  // Find weights from config by metric names
  const { F1, RMSD, κ_SRME } = cps_config

  // Check if any metrics with non-zero weights are missing
  if (
    (F1.weight > 0 && (f1 === undefined || isNaN(f1))) ||
    (RMSD.weight > 0 && (rmsd === undefined || isNaN(rmsd))) ||
    (κ_SRME.weight > 0 && (kappa === undefined || isNaN(kappa)))
  ) return null

  // Skip the calculation if all weights are zero
  const total_weight = F1.weight + RMSD.weight + κ_SRME.weight
  if (total_weight === 0) return 0

  // Calculate weighted sum
  let weighted_sum = 0

  // Add F1 contribution if available and weighted
  if (f1 !== undefined && !isNaN(f1) && F1.weight > 0) {
    weighted_sum += normalize_f1(f1) * F1.weight
  }

  // Add RMSD contribution if available and weighted
  if (rmsd !== undefined && !isNaN(rmsd) && RMSD.weight > 0) {
    weighted_sum += normalize_rmsd(rmsd) * RMSD.weight
  }

  // Add kappa contribution if available and weighted
  if (kappa !== undefined && !isNaN(kappa) && κ_SRME.weight > 0) {
    weighted_sum += normalize_kappa_srme(kappa) * κ_SRME.weight
  }

  // Return weighted average
  return weighted_sum / total_weight
}
