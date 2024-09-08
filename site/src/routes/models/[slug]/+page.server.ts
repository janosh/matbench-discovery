import type { ModelData } from '$lib'
import model_stats from '$lib/model-stats-uniq-protos.json'
import model_schema from '$root/tests/model-schema.yml'
import fs from 'fs'
import { compile as json_to_ts } from 'json-schema-to-typescript'
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

// keep model-schema.d.ts in sync with model-schema.yml (source of truth)
// i.e. use json-schema-to-typescript to auto-convert YAML schema to TypeScript interface
const model_md_type = await json_to_ts(model_schema, `ModelMetadata`, {
  style: { singleQuote: true, semi: false, printWidth: 100 },
})
const dts_out_file = `src/lib/model-schema.d.ts`
fs.writeFileSync(dts_out_file, model_md_type)
