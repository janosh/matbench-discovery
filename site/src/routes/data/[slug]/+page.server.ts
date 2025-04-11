import { DATASETS } from '$lib'
import { error } from '@sveltejs/kit'
import type { PageServerLoad } from './$types'

export const load: PageServerLoad = ({ params }) => {
  const slug = params.slug

  // Find dataset by matching slugs
  const dataset = Object.values(DATASETS).find((dataset) => dataset.slug === slug)

  if (!dataset) {
    throw error(404, { message: `Dataset "${slug}" not found` })
  }

  return { dataset }
}
