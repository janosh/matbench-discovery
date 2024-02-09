import type { ModelData } from '$lib'
import model_stats from '$lib/model-stats.json'
import { compile } from 'mdsvex'
import { dirname } from 'path'

export const load = async () => {
  const files = import.meta.glob(`$root/models/**/*.yml`, {
    eager: true,
    import: `default`,
  }) as Record<string, ModelData>

  // merge computed and static model metadata
  const models = Object.entries(files)
    .filter(([key]) => key.match(/models\/(.+)\/\1.*\.yml/))
    .map(([key, metadata]) => {
      const { model_name } = metadata
      const stats = model_stats[model_name] as ModelData
      if (!stats) console.error(`Missing stats for ${model_name}`)
      return { ...metadata, ...(stats ?? {}), dirname: dirname(key) }
    }) as ModelData[]

  // parse markdown notes to html with mdsvex
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
