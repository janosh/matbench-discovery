import type { ModelData } from '$lib'
import { compile } from 'mdsvex'
import { dirname } from 'path'
import type { PageServerLoad } from './$types'
import model_stats from './2023-02-03-model-stats.json'

export const load: PageServerLoad = async () => {
  const yml = import.meta.glob(`$root/models/**/metadata.yml`, {
    eager: true,
  })

  // merge computed and static model metadata
  const models: ModelData[] = Object.entries(yml).flatMap(([key, module]) => {
    let metadata = module.default as ModelData[]

    if (!Array.isArray(metadata)) metadata = [metadata]

    return metadata.map((md) => {
      const computed = model_stats[md.model_name]
      if (!computed) console.error(`Missing stats for ${md.model_name}`)
      return { ...md, ...(computed ?? {}), dir: dirname(key) }
    })
  })

  // markdown notes to html
  for (const { model_name, notes } of models) {
    if (!notes) continue
    for (const [key, note] of Object.entries(notes)) {
      const out = await compile(note)
      if (!out?.code) {
        console.error(`Failed to compile model note ${model_name}/${key}`)
      } else notes[key] = out.code
    }
  }

  return { models }
}
