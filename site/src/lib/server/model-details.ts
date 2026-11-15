// Model YAML fields the client bundles omit (ServerOnlyModelKey). Server-only, so SvelteKit
// rejects any client import and the full YAML (kept by the SSR transform) stays here.
import type { ModelMetadata } from '#lib/schema/model.js'
import type { ModelDetails } from '#lib/types.js'

// whole modules, not `import: 'default'`: the production build compiled that option to
// the module namespaces (so every `.model_key` read undefined and prerender failed)
const MODEL_YAMLS = Object.values(
  import.meta.glob<{ default: ModelMetadata }>(`$root/models/[^_]**/[^_]*.yml`, {
    eager: true,
  }),
).map((module) => module.default)

export function model_details(model_key: string): ModelDetails {
  const model = MODEL_YAMLS.find((meta) => meta.model_key === model_key)
  if (!model) throw new Error(`No model YAML for ${model_key}`)
  const { notes, environment, trained_by } = model
  return { notes, environment, trained_by }
}
