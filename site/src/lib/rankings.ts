// Per-metric leaderboard ranks for a single model, shown as a "report card" on
// model detail pages. Values are read via the same label paths the metrics table
// uses, so ranks always agree with the leaderboard.
import { ALL_METRICS, MD_METRICS } from '$lib/labels'
import { get_nested_number, is_finite_num, label_data_path } from '$lib/metrics'
import type { Label, ModelData } from '$lib/types'

export interface MetricRank {
  metric: Label
  rank: number // 1-based competition ranking (ties share the better rank)
  n_models: number // number of models with a finite value for this metric
  value: number
}

// headline metrics to rank models by, with the page where each leaderboard lives.
// labels carry their task name for context on model detail pages
export const RANKED_METRICS: (Label & { rank_href: string })[] = [
  { ...ALL_METRICS.CPS, rank_href: `/` },
  { ...ALL_METRICS.F1, label: `Discovery F1`, rank_href: `/` },
  { ...ALL_METRICS.MAE, label: `Discovery MAE`, rank_href: `/` },
  { ...ALL_METRICS.RMSD, label: `Geo Opt RMSD`, rank_href: `/tasks/geo-opt` },
  {
    ...ALL_METRICS.κ_SRME,
    label: `Phonons κ<sub>SRME</sub>`,
    rank_href: `/tasks/phonons`,
  },
  { ...MD_METRICS.md_combined_score, label: `MD CMDS`, rank_href: `/tasks/md` },
]

const metric_value = (model: ModelData, metric: Label): number | undefined => {
  const value = get_nested_number(model, label_data_path(metric))
  return is_finite_num(value) ? value : undefined
}

// Compute the model's rank among `models` for each metric it has a value for.
// The model is resolved from `models` by key so computed-on-the-fly scores
// (CPS/CMDS, which only exist on the shared MODELS entries) are picked up even
// when the caller holds a separately-loaded copy of the model metadata.
export function model_metric_ranks<M extends Label>(
  model_key: string,
  models: ModelData[],
  metrics: readonly M[],
): (MetricRank & { metric: M })[] {
  const model = models.find((entry) => entry.model_key === model_key)
  if (!model) return []

  return metrics.flatMap((metric) => {
    const value = metric_value(model, metric)
    if (value === undefined) return []
    const all_values = models
      .map((entry) => metric_value(entry, metric))
      .filter(is_finite_num)
    const n_better = all_values.filter((other) =>
      metric.better === `lower` ? other < value : other > value,
    ).length
    return [{ metric, rank: n_better + 1, n_models: all_values.length, value }]
  })
}
