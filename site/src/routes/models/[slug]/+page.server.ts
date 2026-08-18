import { MODELS } from '$lib'
import { read_md_per_system } from '$lib/server/md'
import { error, redirect } from '@sveltejs/kit'
import type { EntryGenerator, PageServerLoad } from './$types'

export const entries: EntryGenerator = () =>
  MODELS.flatMap(({ model_key, model_key_aliases = [] }) => [
    { slug: model_key },
    ...model_key_aliases.map((slug) => ({ slug })),
  ])

export const load: PageServerLoad = async ({ params, url }) => {
  const model = MODELS.find(
    ({ model_key, model_key_aliases = [] }) =>
      model_key === params.slug || model_key_aliases.includes(params.slug),
  )

  if (!model) {
    error(404, { message: `Model "${params.slug}" not found` })
  }
  if (model.model_key !== params.slug) {
    redirect(308, `/models/${model.model_key}${url.search}`)
  }

  return { model, md_per_system: await read_md_per_system(model) }
}
