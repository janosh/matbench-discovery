import type { ModelData } from '$lib'
import model_stats from '$lib/model-stats-uniq-protos.json'
import model_schema from '$root/tests/model-schema.yml'
import { prettier as prettier_config } from '$site/package.json'
import fs from 'fs'
import { compile as json_to_ts } from 'json-schema-to-typescript'
import { compile } from 'mdsvex'
import { dirname } from 'path'

export const load = async () => {
  const files = import.meta.glob(`$root/models/[^_]**/[^_]*.yml`, {
    eager: true,
    import: `default`,
  }) as Record<string, ModelData>

  // merge computed and static model metadata
  const models = Object.entries(files)
    .filter(
      // ignore models that aren't completed
      ([_key, metadata]) => (metadata?.status ?? `complete`) == `complete`,
    )
    .map(([key, metadata]) => {
      const { model_name } = metadata
      const stats = model_stats[model_name] as ModelData
      if (!stats) console.trace(`Missing stats for ${model_name}`)
      return { ...metadata, ...(stats ?? {}), dirname: dirname(key) }
    }) as ModelData[]

  // parse markdown notes to html with mdsvex
  for (const { model_name, notes } of models) {
    if (!notes) continue
    for (const [key, note] of Object.entries(notes)) {
      const out = await compile(note)
      if (!out?.code) {
        console.trace(`Failed to compile model note ${model_name}/${key}`)
        // remove outer p tags
      } else notes[key] = out.code.replace(/^\s<p>(.*)<\/p>\s$/, `$1`)
    }
  }

  return { models }
}

// keep model-schema.d.ts in sync with model-schema.yml (source of truth)
// i.e. use json-schema-to-typescript to auto-convert YAML schema to TypeScript interface
const model_metadata_ts = await json_to_ts(model_schema, `ModelMetadata`, {
  style: prettier_config,
})
// prettier format model_md_type
const dts_out_file = `src/lib/model-schema.d.ts`
fs.writeFileSync(dts_out_file, model_metadata_ts)
