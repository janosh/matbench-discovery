// Per-metric leaderboard ranks for a single model, shown as a "report card" on
// model detail pages. Values are read via the same label paths the metrics table
// uses, so ranks always agree with the leaderboard.
import { ALL_METRICS, MD_METRICS } from '$lib/labels'
import { get_nested_number, is_finite_num, label_data_path } from '$lib/metrics'
import type { Label, ModelData } from '$lib/types'

// headline metrics to rank models by, with the page where each leaderboard lives.
// labels carry their task name for context on model detail pages
export const RANKED_METRICS: (Label & { rank_href: string })[] = [
  { ...ALL_METRICS.CPS, rank_href: `/` },
  { ...ALL_METRICS.F1, label: `Discovery F1`, rank_href: `/tasks/discovery` },
  { ...ALL_METRICS.RMSD, label: `Geo Opt RMSD`, rank_href: `/tasks/geo-opt` },
  {
    ...ALL_METRICS.κ_SRME,
    label: `Phonons κ<sub>SRME</sub>`,
    rank_href: `/tasks/phonons`,
  },
  { ...MD_METRICS.md_combined_score, label: `MD CMDS`, rank_href: `/tasks/md` },
  {
    ...ALL_METRICS.diatomics_combined_score,
    label: `Diatomics CDS`,
    rank_href: `/tasks/diatomics`,
  },
]

// Compute the model's rank among `models` for each metric it has a value for.
// Resolve the model from the supplied cohort so computed-on-the-fly scores are used.
export function model_metric_ranks<M extends Label>(
  model_key: string,
  models: ModelData[],
  metrics: readonly M[],
) {
  const model = models.find((entry) => entry.model_key === model_key)
  if (!model) return []

  return metrics.flatMap((metric) => {
    const path = label_data_path(metric)
    const value = get_nested_number(model, path)
    if (!is_finite_num(value)) return []
    const all_values = models
      .map((entry) => get_nested_number(entry, path))
      .filter(is_finite_num)
    const n_better = all_values.filter((other) =>
      metric.better === `lower` ? other < value : other > value,
    ).length
    // rank is 1-based competition ranking (ties share the better rank);
    // n_models counts models with a finite value for this metric
    return [{ metric, rank: n_better + 1, n_models: all_values.length, value }]
  })
}
