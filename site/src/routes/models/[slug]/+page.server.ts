import type { ModelData } from '$lib'
import model_stats from '$lib/model-stats-uniq-protos.json'
import { compile } from 'mdsvex'
import { dirname } from 'path'

export const load = async ({ params }) => {
  const files = import.meta.glob(`$root/models/**/*.yml`, {
    eager: true,
    import: `default`,
  }) as Record<string, ModelData>

  // merge performance metrics and static model metadata
  const [path, metadata] = Object.entries(files).filter(([key]) =>
    key.endsWith(`${params.slug}.yml`),
  )?.[0] as [string, ModelData]

  metadata.dirname = dirname(path)
  const { model_name } = metadata

  // parse markdown notes to html with mdsvex
  if (metadata.notes) {
    for (const [key, note] of Object.entries(metadata.notes)) {
      const out = await compile(note)
      if (!out?.code) {
        console.trace(`Failed to compile model note ${model_name}/${key}`)
        // remove outer p tags
      } else metadata.notes[key] = out.code.replace(/^\s<p>(.*)<\/p>\s$/, `$1`)
    }
  }

  const stats = model_stats[model_name] as ModelData
  if (!stats) console.trace(`Missing stats for ${model_name}`)

  return { model: { ...metadata, ...(stats ?? {}) } }
}
