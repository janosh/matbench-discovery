import type { ModelData, ModelMetadata } from '$lib'
import { compile } from 'mdsvex'
import { dirname } from 'path'
import type { PageServerLoad } from './$types'
import model_stats from './2023-01-23-model-stats.json'

export const load: PageServerLoad = async () => {
  const yml = import.meta.glob(`$root/models/**/metadata.yml`, {
    eager: true,
  })

  // merge computed and static model metadata
  const models: [string, ModelData][] = Object.entries(yml).map(
    ([key, module]) => {
      const metadata = module.default as ModelMetadata
      const computed = model_stats[metadata.model_name] ?? {}
      return [dirname(key), { ...metadata, ...computed }]
    }
  )

  // markdown notes to html
  for (const [name, { notes }] of models) {
    if (!notes) continue
    for (const [key, note] of Object.entries(notes)) {
      const out = await compile(note)
      if (!out?.code) {
        console.error(`Failed to compile model note ${name}/${key}`)
      } else notes[key] = out.code
    }
  }

  return { models }
}
