import rehypeStringify from 'rehype-stringify'
import remarkParse from 'remark-parse'
import remarkRehype from 'remark-rehype'
import { unified } from 'unified'
import type { ModelMetadata } from './model-schema'
import type { ModelData } from './types'

export { default as CaptionedMetricsTable } from './CaptionedMetricsTable.svelte'
export { default as Footer } from './Footer.svelte'
export { default as AuthorBrief } from './ModelAuthor.svelte'
export { default as ModelCard } from './ModelCard.svelte'
export { default as Nav } from './Nav.svelte'
export { default as PtableHeatmap } from './PtableHeatmap.svelte'
export { default as PtableInset } from './PtableInset.svelte'
export { default as References } from './References.svelte'
export * from './types'

export function model_is_compliant(metadata: ModelMetadata): boolean {
  if ((metadata.openness ?? `OSOD`) != `OSOD`) return false

  const allowed_sets = [`MP 2022`, `MPtrj`, `MPF`, `MP Graphs`]

  return metadata.training_set.every((itm) => allowed_sets.includes(itm))
}

const model_metadata_files = import.meta.glob(`$root/models/[^_]**/[^_]*.yml`, {
  eager: true,
  import: `default`,
}) as Record<string, ModelData>

// merge computed and static model metadata
export const MODEL_METADATA = Object.entries(model_metadata_files)
  .filter(
    // ignore models that aren't completed
    ([_key, metadata]) => (metadata?.status ?? `complete`) == `complete`,
  )
  .map(([key, metadata]) => {
    return {
      ...metadata,
      dirname: key.split(`/`)[2],
      metadata_file: key.replace(/^../, ``),
    }
  }) as ModelData[]

// parse markdown notes to html with remark/rehype
for (const { model_name, notes } of MODEL_METADATA) {
  if (!notes) continue
  for (const [key, note] of Object.entries(notes)) {
    const out = unified()
      .use(remarkParse)
      .use(remarkRehype)
      .use(rehypeStringify)
      .processSync(note as string)
    if (!out) {
      console.trace(`Failed to compile model note ${model_name}/${key}`)
      // remove outer p tags
    } else notes[key] = out
  }
}
