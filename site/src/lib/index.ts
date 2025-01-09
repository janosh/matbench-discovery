import { default as data_files } from '$pkg/data-files.yml'
import modeling_tasks from '$pkg/modeling-tasks.yml'
import rehypeStringify from 'rehype-stringify'
import remarkParse from 'remark-parse'
import remarkRehype from 'remark-rehype'
import { unified } from 'unified'
import type { ModelMetadata } from './model-schema'
import type { ModelData } from './types'

export { default as TRAINING_SETS } from '$data/training-sets.yml'
export { default as Footer } from './Footer.svelte'
export { default as GeoOptMetricsTable } from './GeoOptMetricsTable.svelte'
export { default as HeatmapTable } from './HeatmapTable.svelte'
export { default as MetricsTable } from './MetricsTable.svelte'
export { default as AuthorBrief } from './ModelAuthor.svelte'
export { default as ModelCard } from './ModelCard.svelte'
export { default as Nav } from './Nav.svelte'
export { default as PtableHeatmap } from './PtableHeatmap.svelte'
export { default as PtableInset } from './PtableInset.svelte'
export { default as References } from './References.svelte'
export * from './types'
export { data_files }

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
      metadata_file: key.replace(/^..\//, ``),
    }
  }) as ModelData[]

const md_parser = unified().use(remarkParse).use(remarkRehype).use(rehypeStringify)
function md_to_html(md: string): Uint8Array | string {
  return md_parser.processSync(md)?.value
}

// parse markdown notes to html with remark/rehype
for (const { notes, metadata_file } of MODEL_METADATA) {
  if (!notes) continue
  if (!notes.html) notes.html = {}

  for (const [key, note] of Object.entries(notes)) {
    // skip if note was already parsed to HTML or is not a string
    if (typeof note !== `string` || key in notes.html) continue

    const html_note = md_to_html(note)

    if (html_note) notes.html[key] = html_note
    else console.error(`${metadata_file}: Failed to compile note '${key}'\n`)
  }
}

for (const key of Object.keys(data_files).filter((key) => !key.startsWith(`_`))) {
  data_files[key].html = md_to_html(
    data_files[key].description + `\n\n${data_files._links}`,
  )
}

export const all_higher_better = Object.values(modeling_tasks).flatMap(
  (model_task) => model_task.metrics.higher_is_better,
)
export const all_lower_better = Object.values(modeling_tasks).flatMap(
  (model_task) => model_task.metrics.lower_is_better,
)
export function get_metric_rank_order(metric: string) {
  if (all_higher_better.includes(metric)) return `higher`
  if (all_lower_better.includes(metric)) return `lower`
  return null
}
