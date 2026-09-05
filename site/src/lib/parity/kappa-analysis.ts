import type kappa_analysis_data from '$figs/kappa-103-analysis.jsonl'

export type KappaAnalysis = typeof kappa_analysis_data

let analysis_promise: Promise<KappaAnalysis> | undefined

// Lazily load and cache the shared κ-103 per-material analysis payload.
export const load_kappa_analysis = (): Promise<KappaAnalysis> =>
  (analysis_promise ??= import(`$figs/kappa-103-analysis.jsonl`).then(
    (module) => module.default,
  ))

// Map material IDs to one model's per-material κ_SRME values.
export async function load_kappa_srme_map(
  model_key: string,
): Promise<Map<string, number | null> | undefined> {
  const analysis = await load_kappa_analysis()
  const model_analysis = analysis.models.find(
    (model_data) => model_data.model_key === model_key,
  )
  if (!model_analysis) return undefined
  return new Map(
    analysis.material_ids.map((material_id, idx) => [
      material_id,
      model_analysis.srme[idx],
    ]),
  )
}
