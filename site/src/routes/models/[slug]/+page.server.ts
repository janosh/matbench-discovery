import type { ModelData } from '$lib'
import model_stats from '$lib/model-stats-uniq-protos.json'
import { compile } from 'mdsvex'
import { dirname } from 'path'

export const load = async ({ params }) => {
  const files = import.meta.glob(`$root/models/**/*.yml`, {
    eager: true,
    import: `default`,
  }) as Record<string, ModelData>

  const all_metadata = Object.fromEntries(
    Object.entries(files).map(([path, metadata]) => [
      metadata.model_name.replaceAll(` `, `-`).toLowerCase(),
      { ...metadata, dirname: dirname(path) },
    ]),
  )

  // merge performance metrics and static model metadata
  const metadata = all_metadata[params.slug]

  if (!metadata) {
    throw `Failed to load model metadata for ${params.slug}, available: ${Object.keys(all_metadata)}`
  }

  const { model_name } = metadata

  // parse markdown notes to html with mdsvex
  if (metadata.notes) {
    for (const [key, note] of Object.entries(metadata.notes)) {
      if (typeof note !== `string`) continue
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
