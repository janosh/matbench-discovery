import { default as DATASETS } from '$data/datasets.yml'
import { default as data_files } from '$pkg/data-files.yml'
import rehypeStringify from 'rehype-stringify'
import remarkParse from 'remark-parse'
import remarkRehype from 'remark-rehype'
import { unified } from 'unified'
import { MODELS } from './models.svelte'

export { default as Footer } from './Footer.svelte'
export { default as GeoOptMetricsTable } from './GeoOptMetricsTable.svelte'
export { default as Logo } from './Logo.svelte'
export { default as MetricsTable } from './MetricsTable.svelte'
export { default as AuthorBrief } from './ModelAuthor.svelte'
export { default as ModelCard } from './ModelCard.svelte'
export { MODELS } from './models.svelte'
export { default as OrgLogos } from './OrgLogos.svelte'
export { default as PtableHeatmap } from './PtableHeatmap.svelte'
export { default as PtableInset } from './PtableInset.svelte'
export { default as References } from './References.svelte'
export { default as SelectToggle } from './SelectToggle.svelte'
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
export const md_to_html = (md: string): string =>
  String(md_parser.processSync(md)?.value ?? ``)

// Function to slugify text for URLs
export const slugify = (text: string): string =>
  text.toLowerCase().replaceAll(/[\s_]+/g, `-`)

// Convert array types to strings and handle missing values
export function arr_to_str(value: unknown): string {
  if (value === null || value === undefined || value === ``) return `n/a`
  if (Array.isArray(value)) return value.join(`, `)
  return JSON.stringify(value)
}

// Process datasets to add slugs and convert descriptions to HTML
for (const [key, dataset] of Object.entries(DATASETS)) {
  dataset.slug = slugify(key)
  dataset.description_html = md_to_html(dataset.description)
}

// Parse markdown notes to html with remark/rehype
for (const { notes, metadata_file } of MODELS) {
  if (!notes) continue
  notes.html ??= {}

  for (const [key, note] of Object.entries(notes)) {
    // Skip if note was already parsed to HTML or is not a string
    if (typeof note !== `string` || key in notes.html) continue

    const html_note = md_to_html(note)

    if (html_note) notes.html[key] = html_note
    else console.error(`${metadata_file}: Failed to compile note '${key}'\n`)
  }
}

// oxlint-disable-next-line typescript/dot-notation -- `_links` is an external data-files.yml field; dot access trips no-underscore-dangle
const data_file_links = data_files[`_links`]
if (typeof data_file_links !== `string`) {
  throw new TypeError(`data-files.yml: _links must be a string`)
}

for (const [key, entry] of Object.entries(data_files)) {
  if (key.startsWith(`_`) || typeof entry !== `object`) continue
  entry.html = md_to_html(`${entry.description}\n\n${data_file_links}`)
}

// Format date string into human-readable format
export const format_date = (
  date: string | number,
  options?: Intl.DateTimeFormatOptions,
): string =>
  new Date(date).toLocaleDateString(undefined, {
    year: `numeric`,
    month: `short`,
    day: `numeric`,
    ...options,
  })

// Compare two models by submission date, newest first (for Array.toSorted)
export const by_date_added_desc = (
  model_1: { date_added: string },
  model_2: { date_added: string },
): number =>
  new Date(model_2.date_added).getTime() - new Date(model_1.date_added).getTime()
