import type { ModelData, ModelMetadata } from '$lib/types'
import { dirname } from 'path'
import type { PageServerLoad } from './$types'
import model_stats from './2023-01-23-model-stats.json'

export const load: PageServerLoad = async () => {
  const yml = import.meta.glob(`$root/models/**/metadata.yml`, {
    eager: true,
  })
  const models: [string, ModelData][] = Object.entries(yml).map(
    ([key, module]) => {
      const metadata = module.default as ModelMetadata
      const computed = model_stats[metadata.model_name] ?? {}
      return [dirname(key), { ...metadata, ...computed }]
    }
  )

  return { models }
}
