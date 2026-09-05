import { strip_html } from 'matterviz'
import { csv_line } from 'svelte-widgets/csv'
import { ALL_METRICS, HYPERPARAMS, METADATA_COLS } from './labels'

export const heatmap_class = `heatmap`

type ExportOptions = { discovery_set?: string }
type ExportResult = { filename: string; url: string }

// Headers and column indices to export, excluding SVG icon columns (Org and Links)
// and the structural rank (#) column (just 1..N under the current sort)
function get_export_columns(table_el: Element): { headers: string[]; indices: number[] } {
  const header_rows = table_el.querySelectorAll(`thead tr`)
  const main_header_row = header_rows.item(header_rows.length - 1)
  const headers: string[] = []
  const indices: number[] = []
  const all_headers = [...(main_header_row?.querySelectorAll(`th`) ?? [])]
  all_headers.forEach((th, col_idx) => {
    if (th.classList.contains(`row-num-col`)) return
    const header_text = th.textContent?.replaceAll(/[↑↓]/g, ``).trim() || ``
    if (header_text !== `Org` && header_text !== `Links`) {
      headers.push(header_text)
      indices.push(col_idx)
    }
  })
  return { headers, indices }
}

function extract_table_data(): { headers: string[]; rows: (string | number)[][] } {
  const table_el = document.querySelector(`.${heatmap_class}`)

  if (!table_el) {
    throw new Error(`Table element not found for export`)
  }

  const { headers, indices } = get_export_columns(table_el)

  const rows = [...table_el.querySelectorAll(`tbody tr`)].map((row) => {
    const all_cells = [...row.querySelectorAll(`td`)]
    const row_data: (string | number)[] = []

    indices.forEach((col_index) => {
      const cell = all_cells[col_index]
      if (cell) row_data.push(format_cell(cell, headers[row_data.length]))
    })

    return row_data
  })
  return { headers, rows }
}

// Convert a table cell to its exported value: numbers (from data-sort-value or
// text) are formatted per the column's label spec; otherwise return cleaned text
function format_cell(cell: Element, header: string): string | number {
  // Prefer data-sort-value attribute (holds the unformatted number). HeatmapTable leaves
  // it null on the <td> for HTML-string cells, whose raw value metrics.ts puts on an
  // inner <span data-sort-value> instead (model_params, training set, r_cut, ...)
  const sort_value =
    cell.getAttribute(`data-sort-value`) ??
    cell.querySelector(`[data-sort-value]`)?.getAttribute(`data-sort-value`)
  if (sort_value && sort_value !== `null`) {
    const num_value = Number(sort_value)
    return isNaN(num_value) ? sort_value : format_value_for_export(num_value, header)
  }

  // Fall back to text content, parsing as a number when it looks like one
  const text_content = cell.textContent?.trim() || ``
  const num_value = Number(text_content)
  if (!isNaN(num_value) && text_content !== `` && text_content !== `n/a`) {
    return format_value_for_export(num_value, header)
  }
  return text_content
    .replaceAll(`&lt;`, `<`)
    .replaceAll(`&gt;`, `>`)
    .replaceAll(`&amp;`, `&`)
    .replaceAll(/\s+/g, ` `)
    .trim()
}

function format_number(value: number, format?: string): number | string {
  // SI/scientific formats use exponential notation for large numbers
  if ((format?.includes(`~s`) || format?.includes(`.3s`)) && Math.abs(value) >= 1000)
    return value.toExponential(2)

  // Parse decimal places from format strings like '.2f' (default 3)
  const n_decimals = Number(format?.match(/\.(?<decimals>\d)f/)?.groups?.decimals ?? 3)
  const factor = 10 ** n_decimals
  return Math.round(value * factor) / factor
}

// Format a value per the format spec of the label whose text matches the header
function format_value_for_export(value: number, header: string): number | string {
  const all_labels = { ...ALL_METRICS, ...METADATA_COLS, ...HYPERPARAMS }
  const clean_header = strip_html(header).trim()

  const label = Object.values(all_labels).find(
    (lbl) =>
      (lbl.label && strip_html(lbl.label).trim() === clean_header) ||
      (lbl.key && strip_html(lbl.key).trim() === clean_header),
  )

  return format_number(value, label?.format)
}

export function generate_csv({
  discovery_set = `unique_prototypes`,
}: ExportOptions): ExportResult | null {
  try {
    const { headers, rows } = extract_table_data()

    const csv_content = [headers, ...rows].map(csv_line).join(`\n`)

    const blob = new Blob([csv_content], { type: `text/csv;charset=utf-8;` })
    const date = new Date().toISOString().split(`T`)[0]
    const discovery = discovery_set.replaceAll(`_`, `-`)
    const filename = `matbench-${discovery}-${rows.length}models-${date}.csv`
    const url = URL.createObjectURL(blob)
    const anchor = document.createElement(`a`)
    anchor.href = url
    anchor.download = filename
    anchor.click()
    setTimeout(() => URL.revokeObjectURL(url), 100)
    return { filename, url }
  } catch (error) {
    console.error(`Error generating CSV:`, error)
    return null
  }
}
