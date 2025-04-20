import type { CombinedMetricConfig } from './types'

export const RMSD_BASELINE = 0.15 // baseline for poor performance given worst performing model at time of writing is M3GNet at 0.1117

export const DEFAULT_CPS_CONFIG: CombinedMetricConfig = {
  label: `CPS`,
  name: `Combined Performance Score`,
  key: `cps`,
  description: `Combined Performance Score averages discovery (F1), structure optimization (RMSD), and phonon performance (κ<sub>SRME</sub>) according to user-defined weights`,
  range: [0, 1],
  parts: {
    F1: {
      path: `discovery.unique_prototypes.F1`,
      label: `F1`,
      description: `F1 score for stable/unstable material classification (discovery task)`,
      weight: 0.5,
      range: [0, 1],
      better: `higher`,
    },
    kappa_SRME: {
      path: `phonons.kappa_103.κ_SRME`,
      label: `κ<sub>SRME</sub>`,
      svg_label: `κ<tspan baseline-shift='-0.4em' font-size='0.8em'>SRME</tspan>`,
      description: `Symmetric relative mean error of predicted lattice thermal conductivity`,
      weight: 0.4,
      range: [0, 2],
      better: `lower`,
    },
    RMSD: {
      path: `discovery.unique_prototypes.RMSD`,
      label: `RMSD`,
      description: `Root mean square displacement for crystal structure optimization`,
      weight: 0.1,
      range: [0, RMSD_BASELINE],
      better: `lower`,
    },
  },
} as const

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

// kappa_SRME is symmetric relative mean error, with range [0,2] by definition
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
  config: CombinedMetricConfig, // weights for each metric
): number | null {
  // Find weights from config by metric names
  const { F1, RMSD, kappa_SRME } = config.parts

  // Check if any metrics with non-zero weights are missing
  if (
    (F1.weight > 0 && (f1 === undefined || isNaN(f1))) ||
    (RMSD.weight > 0 && (rmsd === undefined || isNaN(rmsd))) ||
    (kappa_SRME.weight > 0 && (kappa === undefined || isNaN(kappa)))
  ) {
    return null
  }

  // Skip the calculation if all weights are zero
  const total_weight = F1.weight + RMSD.weight + kappa_SRME.weight
  if (total_weight === 0) {
    return 0
  }

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
  if (kappa !== undefined && !isNaN(kappa) && kappa_SRME.weight > 0) {
    weighted_sum += normalize_kappa_srme(kappa) * kappa_SRME.weight
  }

  // Return weighted average
  return weighted_sum / total_weight
}
