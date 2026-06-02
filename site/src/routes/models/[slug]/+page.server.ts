import { MODELS } from '$lib'
import { error, redirect } from '@sveltejs/kit'
import type { EntryGenerator, PageServerLoad } from './$types'

// model_keys were switched from underscores to hyphens. These legacy underscore slugs
// are prerendered as permanent redirects to their kebab-cased pages. 2026-06-02
const LEGACY_SLUGS = [`equiformer_v3_mp`, `equiformer_v3_oam`, `grace_2l_oam_l`]

export const load: PageServerLoad = async ({ params }) => {
  const model = MODELS.find((candidate) => candidate.model_key === params.slug)

  if (!model) {
    // schema forbids underscores in model_key, so any underscore slug is a legacy
    // URL -> forward to its kebab-cased equivalent if that resolves to a model
    const kebab = params.slug.replaceAll(`_`, `-`)
    if (MODELS.some((candidate) => candidate.model_key === kebab)) {
      redirect(308, `/models/${kebab}`)
    }
    error(404, { message: `Model "${params.slug}" not found` })
  }

  return { model }
}

export const entries: EntryGenerator = () => LEGACY_SLUGS.map((slug) => ({ slug }))
