import { DATASETS } from '$lib'
import { error } from '@sveltejs/kit'
import type { PageServerLoad } from './$types'

export const load: PageServerLoad = ({ params }) => {
  // Find dataset by matching slugs
  const dataset = Object.values(DATASETS).find((dataset) => dataset.slug === params.slug)

  if (!dataset) {
    throw error(404, { message: `Dataset "${params.slug}" not found` })
  }

  return { dataset }
}
