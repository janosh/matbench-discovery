import type { CombinedMetricConfig, DiscoverySet, HeatmapColumn } from './types'

export const METADATA_COLS: HeatmapColumn[] = [
  { label: `Model`, sticky: true },
  { label: `Training Set`, tooltip: `Size of and link to model training set` },
  { label: `Params`, tooltip: `Number of trainable model parameters` },
  { label: `Targets`, tooltip: `Target property used to train the model` },
  { label: `Date Added`, tooltip: `Submission date to the leaderboard` },
  {
    label: `Links`,
    tooltip: `Model resources: paper, code repository and submission pull request`,
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
    link: `https://github.com/janosh/matbench-discovery/blob/fd1dda6c/data/wbm/compile_wbm_test_set.py#L632-L705`,
  },
  most_stable_10k: {
    title: `10k Most Stable`,
    tooltip: `Metrics computed on the 10k structures predicted to be most stable (different for each model)`,
  },
}

export const [F1_DEFAULT_WEIGHT, RMSD_DEFAULT_WEIGHT, KAPPA_DEFAULT_WEIGHT] = [
  0.5, 0.1, 0.4,
]

export const DEFAULT_COMBINED_METRIC_CONFIG: CombinedMetricConfig = {
  name: `CPS`,
  description: `Combined Performance Score weighting discovery, structure optimization, and phonon performance`,
  weights: [
    {
      metric: `F1`,
      label: `F1`,
      description: `F1 score for stable/unstable material classification (discovery task)`,
      value: F1_DEFAULT_WEIGHT,
    },
    {
      metric: `kappa_SRME`,
      label: `κ<sub>SRME</sub>`,
      description: `Symmetric relative mean error for thermal conductivity prediction (lower is better)`,
      value: KAPPA_DEFAULT_WEIGHT,
    },
    {
      metric: `RMSD`,
      label: `RMSD`,
      description: `Root mean square displacement for crystal structure optimization`,
      value: RMSD_DEFAULT_WEIGHT,
    },
  ],
}

// F1 score is between 0-1 where higher is better (no normalization needed)
function normalize_f1(value: number | undefined): number {
  if (value === undefined || isNaN(value)) return 0
  return value // Already in [0,1] range
}

// RMSD is lower=better, with current models in the range of ~0.02-0.25 Å
// We invert this so that better performance = higher score
function normalize_rmsd(value: number | undefined): number {
  if (value === undefined || isNaN(value)) return 0

  // Fixed reference points for RMSD (in Å)
  const excellent = 0 // Perfect performance (atoms in exact correct positions)
  const baseline = 0.3 // in Å, a reasonable baseline for poor performance given worst performing model at time of writing is AlphaNet-MPTrj at 0.0227 Å

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

/**
 * Calculate a combined score using normalized metrics weighted by importance factors.
 * This uses fixed normalization reference points to ensure score stability when new models are added.
 *
 * Normalization reference points:
 * - F1: Already in [0,1] range, higher is better
 * - RMSD: 0.0Å (perfect) to 0.25Å (baseline), lower is better
 * - κ_SRME: Range [0,2] linearly mapped to [1,0], lower is better
 *
 * @param f1 F1 score for discovery
 * @param rmsd Root mean square displacement in Å
 * @param kappa Symmetric relative mean error for thermal conductivity
 * @param config Configuration with weights for each metric
 * @returns Combined score between 0-1, or NaN if any weighted metric is missing
 */
export function calculate_combined_score(
  f1: number | undefined,
  rmsd: number | undefined,
  kappa: number | undefined,
  config: CombinedMetricConfig,
): number {
  // Find weights from config by metric names
  const f1_weight =
    config.weights.find((w) => w.metric === `F1`)?.value ?? F1_DEFAULT_WEIGHT
  const rmsd_weight =
    config.weights.find((w) => w.metric === `RMSD`)?.value ?? RMSD_DEFAULT_WEIGHT
  const kappa_weight =
    config.weights.find((w) => w.metric === `kappa_SRME`)?.value ?? KAPPA_DEFAULT_WEIGHT

  // Check if any weighted metric is missing - if so, return NaN
  if (
    (f1_weight > 0 && f1 === undefined) ||
    (rmsd_weight > 0 && rmsd === undefined) ||
    (kappa_weight > 0 && kappa === undefined)
  ) {
    return NaN
  }

  // Get normalized metric values
  const normalized_f1 = normalize_f1(f1)
  const normalized_rmsd = normalize_rmsd(rmsd)
  const normalized_kappa = normalize_kappa_srme(kappa)

  // Get available weights and metrics
  const available_metrics = []
  const available_weights = []

  // Only include metrics that are available
  if (f1 !== undefined) {
    available_metrics.push(normalized_f1)
    available_weights.push(f1_weight)
  }

  if (rmsd !== undefined) {
    available_metrics.push(normalized_rmsd)
    available_weights.push(rmsd_weight)
  }

  if (kappa !== undefined) {
    available_metrics.push(normalized_kappa)
    available_weights.push(kappa_weight)
  }

  // If no metrics are available, return 0
  if (available_metrics.length === 0) return 0

  // Normalize weights to sum to 1 based on available metrics
  const weight_sum = available_weights.reduce((sum, w) => sum + w, 0)
  const normalized_weights =
    weight_sum > 0
      ? available_weights.map((w) => w / weight_sum)
      : available_weights.map(() => 1 / available_weights.length)

  // Calculate weighted average
  let score = 0
  for (let i = 0; i < available_metrics.length; i++) {
    score += available_metrics[i] * normalized_weights[i]
  }

  return score
}
