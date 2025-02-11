import type { HeatmapColumn } from './types'

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
