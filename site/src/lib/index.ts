import { default as DATASETS } from '$data/datasets.yml'
import { default as data_files } from '$pkg/data-files.yml'
import rehypeStringify from 'rehype-stringify'
import remarkParse from 'remark-parse'
import remarkRehype from 'remark-rehype'
import { unified } from 'unified'
import { MODELS } from './models.svelte'

export { default as DiatomicCurve } from './DiatomicCurve.svelte'
export { default as DynamicScatter } from './DynamicScatter.svelte'
export { default as Footer } from './Footer.svelte'
export { default as GeoOptMetricsTable } from './GeoOptMetricsTable.svelte'
export { default as HeatmapTable } from './HeatmapTable.svelte'
export { default as IconList } from './IconList.svelte'
export { default as MetricScatter } from './MetricScatter.svelte'
export { default as MetricsTable } from './MetricsTable.svelte'
export { default as AuthorBrief } from './ModelAuthor.svelte'
export { default as ModelCard } from './ModelCard.svelte'
export { MODELS } from './models.svelte'
export { default as Nav } from './Nav.svelte'
export { default as PtableHeatmap } from './PtableHeatmap.svelte'
export { default as PtableInset } from './PtableInset.svelte'
export { default as RadarChart } from './RadarChart.svelte'
export { default as References } from './References.svelte'
export { default as SelectToggle } from './SelectToggle.svelte'
export { default as TableColumnToggleMenu } from './TableColumnToggleMenu.svelte'
export { default as TableControls } from './TableControls.svelte'
export * from './types'
export { data_files, DATASETS }

// Calculate number of days between provided date string and today
export function calculate_days_ago(date_str: string): string {
  if (!date_str) return ``
  const ms_per_day = 1000 * 60 * 60 * 24
  const now = new Date()
  const then = new Date(date_str)

  return ((now.getTime() - then.getTime()) / ms_per_day).toLocaleString(`en-US`, {
    maximumFractionDigits: 0,
  })
}

const md_parser = unified().use(remarkParse).use(remarkRehype).use(rehypeStringify)
export function md_to_html(md: string) {
  return md_parser.processSync(md)?.value as string
}

// Function to slugify text for URLs
export function slugify(text: string): string {
  return text.toLowerCase().replace(/[\s_]+/g, `-`)
}

// convert array types to strings and handle missing values
export function arr_to_str(value: unknown): string {
  if (!value) return `n/a`
  if (Array.isArray(value)) return value.join(`, `)
  return String(value)
}

// Process datasets to add slugs and convert descriptions to HTML
for (const [key, dataset] of Object.entries(DATASETS)) {
  dataset.slug = slugify(key)
  dataset.description_html = md_to_html(dataset.description)
}

// parse markdown notes to html with remark/rehype
for (const { notes, metadata_file } of MODELS) {
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

// Format date string into human-readable format
export function format_date(date: string, options?: Intl.DateTimeFormatOptions): string {
  return new Date(date).toLocaleDateString(undefined, {
    year: `numeric`,
    month: `short`,
    day: `numeric`,
    ...options,
  })
}
