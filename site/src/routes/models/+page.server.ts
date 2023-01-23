import type { ModelMetadata } from '$lib/types'
import { dirname } from 'path'
import type { PageServerLoad } from './$types'
import analysis from './2023-01-23-pred-analysis.json'

export const load: PageServerLoad = async () => {
  const yml = import.meta.glob(`$root/models/**/metadata.yml`, {
    eager: true,
  })
  const models: [string, ModelMetadata][] = Object.entries(yml).map(
    ([key, module]) => {
      const metadata = module.default as ModelMetadata
      const computed = analysis[metadata.model_name] ?? {}
      return [dirname(key), { ...metadata, ...computed }]
    }
  )

  return { models }
}
