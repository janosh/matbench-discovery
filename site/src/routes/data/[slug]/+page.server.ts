import { DATASETS } from '$lib'
import { error } from '@sveltejs/kit'
import type { PageServerLoad } from './$types'

export const load: PageServerLoad = ({ params }) => {
  const dataset = Object.values(DATASETS).find(({ slug }) => slug === params.slug)

  if (!dataset) {
    error(404, { message: `Dataset "${params.slug}" not found` })
  }

  return { dataset }
}
