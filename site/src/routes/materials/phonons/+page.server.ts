import kappa_analysis from '#figs/kappa-103-analysis.jsonl'
import { kappa_material_summaries } from '#lib/materials.js'
import type { PageServerLoad } from './$types'

export const load: PageServerLoad = () => ({
  materials: kappa_material_summaries(kappa_analysis),
  n_models: kappa_analysis.models.length,
})
