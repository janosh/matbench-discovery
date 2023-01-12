import type { ModelMetadata } from '$lib/types'
import { dirname } from 'path'
import type { PageServerLoad } from './$types'

export const load: PageServerLoad = async () => {
  const models: [string, ModelMetadata][] = Object.entries(
    import.meta.glob(`$root/models/**/metadata.yml`, {
      eager: true,
    })
  ).map(([key, module]) => [dirname(key), module.default])

  return { models }
}
