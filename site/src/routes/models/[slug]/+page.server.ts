import { MODELS } from '#lib/models.svelte.js'
import { read_md_per_system } from '#lib/server/predictions.js'
import { error } from '@sveltejs/kit'
import type { EntryGenerator, PageServerLoad } from './$types'

export const entries: EntryGenerator = () =>
  MODELS.map(({ model_key: slug }) => ({ slug }))

export const load: PageServerLoad = async ({ params }) => {
  const model = MODELS.find(({ model_key }) => model_key === params.slug)

  if (!model) {
    error(404, `Model "${params.slug}" not found`)
  }
  return { model_key: model.model_key, md_per_system: await read_md_per_system(model) }
}
