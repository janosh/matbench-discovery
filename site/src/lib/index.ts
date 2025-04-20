import { default as DATASETS } from '$data/datasets.yml'
import { default as data_files } from '$pkg/data-files.yml'
import modeling_tasks from '$pkg/modeling-tasks.yml'
import rehypeStringify from 'rehype-stringify'
import remarkParse from 'remark-parse'
import remarkRehype from 'remark-rehype'
import { unified } from 'unified'
import { MODELS } from './models.svelte'
import type { ModelData } from './types'

export { default as DiatomicCurve } from './DiatomicCurve.svelte'
export { default as DynamicScatter } from './DynamicScatter.svelte'
export { default as Footer } from './Footer.svelte'
export { default as GeoOptMetricsTable } from './GeoOptMetricsTable.svelte'
export { default as HeatmapTable } from './HeatmapTable.svelte'
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

export function model_is_compliant(model: ModelData): boolean {
  if ((model.openness ?? `OSOD`) != `OSOD`) return false

  const allowed_sets = [`MP 2022`, `MPtrj`, `MPF`, `MP Graphs`]

  return model.training_set.every((itm) => allowed_sets.includes(itm))
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

export function get_pred_file_urls(model: ModelData) {
  // get all pred_file_url from model.metrics
  const files: { name: string; url: string }[] = []

  function find_pred_files(obj: object, parent_key = ``) {
    if (!obj || typeof obj !== `object`) return

    for (const [key, val] of Object.entries(obj)) {
      if (key == `pred_file_url` && val && typeof val === `string`) {
        // Look up the label by traversing the modeling_tasks hierarchy
        const pretty_label = get_label_for_key_path(parent_key)
        files.push({ name: pretty_label, url: val })
      } else if (typeof val === `object`) {
        find_pred_files(val, key)
      }
    }
  }

  // Recursively look up labels in the modeling_tasks object
  function get_label_for_key_path(key_path: string): string {
    if (key_path in modeling_tasks) return modeling_tasks[key_path]?.label || key_path

    // Check if it's a subtask by searching all tasks
    for (const task_value of Object.values(modeling_tasks)) {
      if (task_value?.subtasks?.[key_path]) {
        return task_value.subtasks[key_path].label || key_path
      }
    }

    return key_path // Default to key itself if no label is found
  }

  find_pred_files(model.metrics)
  return files
}
