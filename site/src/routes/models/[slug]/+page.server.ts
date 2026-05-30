import { MODELS } from '$lib'
import { error } from '@sveltejs/kit'
import type { PageServerLoad } from './$types'

export const load: PageServerLoad = async ({ params }) => {
  const model = MODELS.find(
    (candidate_model) => candidate_model.model_key === params.slug,
  )

  if (!model) {
    error(404, { message: `Model "${params.slug}" not found` })
  }

  return { model }
}
