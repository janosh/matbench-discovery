import { MODELS } from '$lib'
import { error } from '@sveltejs/kit'
import per_elem_each_errors from '../per-element-each-errors.json'
import type { PageServerLoad } from './$types'

export const load: PageServerLoad = async ({ params }) => {
  const model = MODELS.find((model) => model.model_key === params.slug)

  if (!model) {
    throw error(404, { message: `Model "${params.slug}" not found` })
  }

  if (!(model.model_name in per_elem_each_errors)) {
    const message = `No per-element energy errors found for ${model.model_name}`
    throw error(404, { message })
  }

  // fail site build if energy parity plots are missing
  for (const which_energy of [`e-form`, `each`]) {
    try {
      await import(`$figs/energy-parity/${which_energy}-parity-${model.model_key}.svelte`)
    } catch (exc) {
      throw error(404, { message: (exc as Error).message })
    }
  }

  return { model }
}
