import { diatomic_models, reference_names } from '#lib/server/diatomics.js'
import type { PageServerLoad } from './$types'

// Model keys only: the page resolves them from MODELS so it sorts by live (re-weighted)
// scores rather than stale serialized copies, and fetches each selected model's curves
// from the prerendered curves/<key>.json endpoint.
export const load: PageServerLoad = () => ({
  diatomic_model_keys: diatomic_models.map(({ model_key }) => model_key),
  reference_names,
})
