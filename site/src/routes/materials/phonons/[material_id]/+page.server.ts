import kappa_analysis from '#figs/kappa-103-analysis.jsonl'
import { kappa_material } from '#lib/materials.js'
import { error } from '@sveltejs/kit'
import type { EntryGenerator, PageServerLoad } from './$types'

export const entries: EntryGenerator = () =>
  kappa_analysis.material_ids.map((material_id) => ({ material_id }))

export const load: PageServerLoad = ({ params }) => {
  if (!kappa_analysis.material_ids.includes(params.material_id)) {
    error(404, `Material "${params.material_id}" is not in the κ-103 test set`)
  }
  return { material: kappa_material(kappa_analysis, params.material_id) }
}
