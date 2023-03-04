import type { ModelData } from '$lib'
import { compile } from 'mdsvex'
import { dirname } from 'path'
import type { PageServerLoad } from './$types'
import model_stats from './model-stats.json'

export const load: PageServerLoad = async () => {
  const yml = import.meta.glob(`$root/models/**/metadata.yml`, {
    eager: true,
  })

  // merge computed and static model metadata
  const models: ModelData[] = Object.entries(yml).flatMap(([key, module]) => {
    let metadata = module.default as ModelData[]

    if (!Array.isArray(metadata)) metadata = [metadata]

    return metadata.flatMap(({ model_name: name, ...rest }) => {
      // handle model name as array
      return (Array.isArray(name) ? name : [name]).map((model_name) => {
        const computed = model_stats[model_name]
        if (!computed) console.error(`Missing stats for ${name}`)
        return { ...rest, model_name, ...(computed ?? {}), dir: dirname(key) }
      })
    })
  })

  // markdown notes to html
  for (const { model_name, notes } of models) {
    if (!notes) continue
    for (const [key, note] of Object.entries(notes)) {
      const out = await compile(note)
      if (!out?.code) {
        console.error(`Failed to compile model note ${model_name}/${key}`)
        // remove outer p tags
      } else notes[key] = out.code.replace(/^\s<p>(.*)<\/p>\s$/, `$1`)
    }
  }

  return { models }
}
