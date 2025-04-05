import type { CombinedMetricConfig, DiscoverySet, HeatmapColumn } from './types'

export const METADATA_COLS: HeatmapColumn[] = [
  { label: `Model`, sticky: true },
  { label: `Training Set`, tooltip: `Size of and link to model training set` },
  { label: `Params`, tooltip: `Number of trainable model parameters` },
  { label: `Targets`, tooltip: `Target property used to train the model` },
  {
    label: `Date Added`,
    tooltip: `Submission date to the leaderboard`,
  },
  {
    label: `Links`,
    tooltip: `Model resources: paper, code repository and submission pull request`,
    sortable: false,
  },
]

export const DISCOVERY_METRICS: HeatmapColumn[] = [
  { label: `F1`, tooltip: `Harmonic mean of precision and recall` },
  { label: `DAF`, tooltip: `Discovery acceleration factor` },
  { label: `Prec`, tooltip: `Precision of classifying thermodynamic stability` },
  { label: `Acc`, tooltip: `Accuracy of classifying thermodynamic stability` },
  {
    label: `TPR`,
    tooltip: `True positive rate of classifying thermodynamic stability`,
  },
  {
    label: `TNR`,
    tooltip: `True negative rate of classifying thermodynamic stability`,
  },
  {
    label: `MAE`,
    tooltip: `Mean absolute error of predicting the convex hull distance`,
    style: `border-left: 1px solid black;`,
  },
  {
    label: `RMSE`,
    tooltip: `Root mean squared error of predicting the convex hull distance`,
  },
  { label: `R<sup>2</sup>`, tooltip: `Coefficient of determination` },
]

export const PHONON_METRICS: HeatmapColumn[] = [
  {
    label: `κ<sub>SRME</sub>`,
    tooltip: `Symmetric relative mean error in predicted phonon mode contributions to thermal conductivity κ`,
    style: `border-left: 1px solid black;`,
  },
]

// Define geometry optimization metrics
export const GEO_OPT_METRICS: HeatmapColumn[] = [
  {
    label: `RMSD`,
    tooltip: `Root mean squared displacement between predicted and reference structures after relaxation`,
    style: `border-left: 1px solid black;`,
  },
  {
    label: `Energy Diff`,
    tooltip: `Mean absolute energy difference between predicted and reference structures`,
  },
  {
    label: `Force RMSE`,
    tooltip: `Root mean squared error of forces in predicted structures relative to reference`,
  },
  {
    label: `Stress RMSE`,
    tooltip: `Root mean squared error of stress in predicted structures relative to reference`,
  },
  {
    label: `Max Force`,
    tooltip: `Maximum force component in predicted structures after relaxation`,
  },
]

// Update ALL_METRICS to include GEO_OPT_METRICS
export const ALL_METRICS: HeatmapColumn[] = [
  ...DISCOVERY_METRICS,
  ...PHONON_METRICS,
  ...GEO_OPT_METRICS.slice(0, 1), // Only include RMSD by default, others can be toggled
]

export const DISCOVERY_SET_LABELS: Record<
  DiscoverySet,
  { title: string; tooltip: string; link?: string }
> = {
  full_test_set: {
    title: `Full Test Set`,
    tooltip: `Metrics computed on the full test set including duplicate structure prototypes`,
  },
  unique_prototypes: {
    title: `Unique Prototypes`,
    tooltip: `Metrics computed only on ~215k unique structure prototypes in WBM determined by matching Aflow-style prototype strings.`,
    link: `https://github.com/janosh/matbench-discovery/blob/37baf7986f848/data/wbm/compile_wbm_test_set.py#L640-L654`,
  },
  most_stable_10k: {
    title: `10k Most Stable`,
    tooltip: `Metrics computed on the 10k structures predicted to be most stable (different for each model)`,
  },
}

export const DEFAULT_CPS_CONFIG: CombinedMetricConfig = {
  label: `CPS`,
  name: `Combined Performance Score`,
  key: `cps`,
  description: `Combined Performance Score weights discovery (F1), structure optimization (RMSD), and phonon performance (κ<sub>SRME</sub>)`,
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
      description: `Symmetric relative mean error for thermal conductivity prediction (lower is better)`,
      weight: 0.4,
      range: [0, 2],
      better: `lower`,
    },
    RMSD: {
      path: `discovery.unique_prototypes.RMSD`,
      label: `RMSD`,
      description: `Root mean square displacement for crystal structure optimization`,
      weight: 0.1,
      range: [0, 0.03],
      better: `lower`,
    },
  },
}

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
  const baseline = 0.03 // baseline for poor performance given worst performing model at time of writing is AlphaNet-MPTrj at 0.0227 Å

  // Linear interpolation between fixed points with clamping
  // Inverse mapping since lower RMSD is better
  if (value <= excellent) return 1.0
  if (value >= baseline) return 0.0
  return (baseline - value) / (baseline - excellent)
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
// - RMSD Root mean square displacement in range 0Å (perfect) to 0.03Å (baseline), lower is better
// - κ_SRME symmetric relative mean error for lattice thermal conductivity,
//    range [0,2] linearly mapped to [1,0], lower is better
export function calculate_combined_score(
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
