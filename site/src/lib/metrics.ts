import type { DiscoverySet, HeatmapColumn } from './types'

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

export const METRICS_COLS: HeatmapColumn[] = [
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
  {
    label: `κ<sub>SRME</sub>`,
    tooltip: `Symmetric relative mean error in predicted phonon mode contributions to thermal conductivity κ`,
    style: `border-left: 1px solid black;`,
  },
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
