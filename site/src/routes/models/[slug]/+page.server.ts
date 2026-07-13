import { MODELS } from '$lib'
import { read_md_per_system } from '$lib/server/md'
import { error } from '@sveltejs/kit'
import type { EntryGenerator, PageServerLoad } from './$types'

export const entries: EntryGenerator = () =>
  MODELS.map(({ model_key }) => ({ slug: model_key }))

export const load: PageServerLoad = async ({ params }) => {
  const model = MODELS.find((candidate) => candidate.model_key === params.slug)

  if (!model) {
    error(404, { message: `Model "${params.slug}" not found` })
  }

  return { model, md_per_system: await read_md_per_system(model) }
}
