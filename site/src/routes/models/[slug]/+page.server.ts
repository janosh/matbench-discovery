import per_elem_each_errors from '$figs/per-element-each-errors.json'
import { MODEL_METADATA } from '$lib'
import { error } from '@sveltejs/kit'
import type { PageServerLoad } from './$types'

export const load: PageServerLoad = ({ params }) => {
  const model = MODEL_METADATA.find((model) => model.model_key === params.slug)

  if (!model) {
    throw error(404, { message: `Model "${params.slug}" not found` })
  }

  if (!(model.model_name in per_elem_each_errors)) {
    const message = `No per-element energy errors found for ${model.model_name}`
    throw error(404, { message })
  }

  return { model }
}
