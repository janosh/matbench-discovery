import { toPng, toSvg } from 'html-to-image'
import { is_dark_mode, strip_html } from 'matterviz'
import { ALL_METRICS, HYPERPARAMS, METADATA_COLS } from './labels'

export const heatmap_class = `heatmap`

type ExportOptions = { discovery_set?: string }
type ExportResult = { filename: string; url: string }

// Simple function to generate descriptive filename with current table state
function generate_filename(
  format: string,
  discovery_set: string = `unique_prototypes`,
): string {
  const date = new Date().toISOString().split(`T`)[0]
  const discovery = discovery_set.replaceAll(`_`, `-`)

  // Get current table state for context
  const table_el = document.querySelector(`.${heatmap_class}`)
  const model_count = table_el?.querySelectorAll(`tbody tr`).length ?? 0

  return `matbench-${discovery}-${model_count}models-${date}.${format.toLowerCase()}`
}

// Create a download link and click it
function trigger_download(url: string, filename: string): void {
  const anchor = document.createElement(`a`)
  anchor.href = url
  anchor.download = filename
  anchor.click()
}

// Create an object URL for a blob, trigger its download, then schedule cleanup
function download_blob(blob: Blob, filename: string): ExportResult {
  const url = URL.createObjectURL(blob)
  trigger_download(url, filename)
  setTimeout(() => URL.revokeObjectURL(url), 100)
  return { filename, url }
}

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

// Helper function to create a filtered table clone excluding SVG icon columns
function create_filtered_table_clone(): HTMLElement {
  const table_el = document.querySelector(`.${heatmap_class}`)

  if (!table_el) throw new Error(`Table element not found for export`)

  // Clone the entire table
  const cloned_node = table_el.cloneNode(true)
  if (!(cloned_node instanceof HTMLTableElement)) {
    throw new Error(`Cloned node is not an HTMLTableElement`)
  }
  const table_clone = cloned_node

  // Remove icon columns from all rows (working backwards to maintain indices)
  const { indices: kept_indices } = get_export_columns(table_clone)
  const header_rows = table_clone.querySelectorAll(`thead tr`)
  const main_header_row = header_rows.item(header_rows.length - 1)
  const n_cols = main_header_row?.querySelectorAll(`th`).length ?? 0
  for (let col_idx = n_cols - 1; col_idx >= 0; col_idx--) {
    if (kept_indices.includes(col_idx)) continue
    for (const row of table_clone.querySelectorAll(`tr`)) {
      row.children[col_idx]?.remove()
    }
  }

  return table_clone
}

// Remove HTML comments recursively (Svelte generates many of these)
function remove_comments(node: Node): void {
  if (node.nodeType === Node.COMMENT_NODE) {
    ;(node as ChildNode).remove()
    return
  }
  for (const child of node.childNodes) remove_comments(child)
}

// Clean up Svelte-specific artifacts and problematic elements
function clean_table_for_export(table_clone: HTMLElement): void {
  // Replace problematic and inline elements with their text content (preserve
  // sub/sup for visual formatting), or remove them entirely if empty
  const replaceable = table_clone.querySelectorAll(
    `svg, img, button, a[href], span, small`,
  )
  replaceable.forEach((el) => {
    if (el.textContent?.trim()) {
      const text_node = document.createTextNode(el.textContent.trim())
      el.parentNode?.replaceChild(text_node, el)
    } else {
      el.remove()
    }
  })

  // Remove HTML comments
  remove_comments(table_clone)

  // Simplify class names and remove data attributes
  const all_elements = [table_clone, ...table_clone.querySelectorAll(`*`)]
  for (const el of all_elements) {
    // Remove Svelte-generated dynamic classes
    if (el.className) {
      const classes = el.className.split(` `).filter((cls) => !cls.startsWith(`s-`))
      el.className = classes.join(` `)
    }

    // Remove Svelte-specific data attributes
    for (const attr of el.attributes) {
      if (
        attr.name.startsWith(`data-`) &&
        ![`data-sort-value`, `data-title`].includes(attr.name)
      ) {
        el.removeAttribute(attr.name)
      }
    }

    // Clean up problematic styling that could cause scrollbars and extra width
    if (!(el instanceof HTMLElement)) continue
    if (el.style) {
      for (const prop of `grid-column grid-row grid-area overflow overflow-x overflow-y
          max-width max-height min-width flex flex-grow flex-basis`.split(/\s+/))
        el.style.removeProperty(prop)

      // Ensure elements don't restrict their size
      el.style.overflow = `visible`
      el.style.whiteSpace = `nowrap`

      // For table cells, apply compact sizing
      if (el.tagName === `TD` || el.tagName === `TH`) {
        el.style.width = `auto`
        el.style.padding = `2px 4px`
        el.style.textAlign = `left`
        el.style.verticalAlign = `middle`
      }
    }
  }
}

// Detect current theme and get appropriate colors
function get_theme_colors(): { background: string; text: string } {
  const root_styles = getComputedStyle(document.documentElement)
  const mode = is_dark_mode() ? `dark` : `light`
  const fallback =
    mode === `dark`
      ? { background: `#061e25`, text: `rgb(208, 208, 208)` }
      : { background: `#fefefe`, text: `#1f2937` }

  return {
    background: root_styles.getPropertyValue(`--${mode}-page-bg`) || fallback.background,
    text: root_styles.getPropertyValue(`--${mode}-text`) || fallback.text,
  }
}

// Create export container with proper styling
function create_export_container(table_clone: HTMLElement): HTMLElement {
  const theme_colors = get_theme_colors()

  const container = document.createElement(`div`)
  Object.assign(container.style, {
    backgroundColor: theme_colors.background,
    color: theme_colors.text,
    display: `inline-block`,
    whiteSpace: `nowrap`,
    position: `relative`,
    width: `fit-content`,
    height: `fit-content`,
    margin: `0`,
    padding: `0`,
    // Ensure no scrollbars appear
    overflow: `visible`,
    overflowX: `visible`,
    overflowY: `visible`,
  })

  const inner = document.createElement(`div`)
  Object.assign(inner.style, {
    padding: `0`,
    margin: `0`,
    boxSizing: `border-box`,
    overflow: `visible`,
    width: `fit-content`,
    height: `fit-content`,
  })
  inner.append(table_clone)
  container.append(inner)

  // Ensure the table itself doesn't create scrollbars or extra spacing
  Object.assign(table_clone.style, {
    overflow: `visible`,
    width: `fit-content`,
    height: `auto`,
    margin: `0`,
    border: `none`,
    tableLayout: `auto`,
    borderCollapse: `collapse`,
  })

  return container
}

// Common filter function for html-to-image
const create_export_filter = (node: Node): boolean => {
  const node_name = node.nodeName

  // Skip external stylesheets
  if (
    node_name === `LINK` &&
    node instanceof Element &&
    node.getAttribute(`rel`) === `stylesheet` &&
    node.getAttribute(`href`)?.startsWith(`http`)
  )
    return false

  // Skip problematic elements
  if ([`SVG`, `IMG`].includes(node_name)) return false
  return true
}

// Log detailed error information
function log_export_error(
  error: unknown,
  format: string,
  context?: Record<string, unknown>,
): void {
  console.error(`Error during ${format} generation:`, {
    error,
    errorType: typeof error,
    errorConstructor: error?.constructor?.name,
    message: error instanceof Error ? error.message : `No message available`,
    stack: error instanceof Error ? error.stack : `No stack available`,
    eventType: error instanceof Event ? error.type : `Not an Event`,
    ...context,
  })
}

// Shared SVG/PNG export: clone + clean table, render off-screen, encode, download
async function generate_image(
  format: `svg` | `png`,
  { discovery_set = `unique_prototypes` }: ExportOptions,
): Promise<ExportResult | null> {
  let container: HTMLElement | undefined
  let table_clone: HTMLElement | undefined
  try {
    table_clone = create_filtered_table_clone()
    clean_table_for_export(table_clone)

    container = create_export_container(table_clone)
    if (format === `png`) {
      container.style.fontFamily = `Arial, sans-serif`
      container.style.fontSize = `14px`
    }
    document.body.append(container)

    const filename = generate_filename(format, discovery_set)

    if (format === `svg`) {
      const svg_data_url = await toSvg(container, {
        backgroundColor: container.style.backgroundColor,
        skipFonts: true,
        filter: create_export_filter,
      })
      const blob = await fetch(svg_data_url).then((res) => res.blob())
      return download_blob(blob, filename)
    }

    // Generate PNG with precise dimensions
    const container_rect = container.getBoundingClientRect()
    const png_data_url = await toPng(container, {
      backgroundColor: container.style.backgroundColor,
      pixelRatio: 2, // High resolution
      width: Math.ceil(container_rect.width),
      height: Math.ceil(container_rect.height),
      skipFonts: true,
      quality: 0.95,
      filter: create_export_filter,
      style: { overflow: `visible`, margin: `0`, padding: `0` },
      imagePlaceholder: `data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQDwAEhQGAhKmMIQAAAABJRU5ErkJggg==`,
    })
    trigger_download(png_data_url, filename)
    return { filename, url: png_data_url }
  } catch (error) {
    const context =
      container && table_clone
        ? {
            containerHTML: `${container.outerHTML.slice(0, 1000)}...`,
            tableCloneHTML: `${table_clone.outerHTML.slice(0, 1000)}...`,
          }
        : undefined
    log_export_error(error, format.toUpperCase(), context)
    return null
  } finally {
    container?.remove()
  }
}

// Function to generate SVG from the metrics table
export const generate_svg = (options: ExportOptions): Promise<ExportResult | null> =>
  generate_image(`svg`, options)

// Function to generate PNG from the metrics table
export const generate_png = (options: ExportOptions): Promise<ExportResult | null> =>
  generate_image(`png`, options)

// Helper function to extract table data excluding SVG icon columns
function extract_table_data(): { headers: string[]; rows: (string | number)[][] } {
  const table_el = document.querySelector(`.${heatmap_class}`)

  if (!table_el) {
    throw new Error(`Table element not found for export`)
  }

  const { headers, indices } = get_export_columns(table_el)

  // Extract data rows
  const data_rows = table_el.querySelectorAll(`tbody tr`)
  const formatted_rows: (string | number)[][] = []

  data_rows.forEach((row) => {
    const all_cells = [...row.querySelectorAll(`td`)]
    const row_data: (string | number)[] = []

    // Only include cells from columns we want to export
    indices.forEach((col_index) => {
      const cell = all_cells[col_index]
      if (cell) row_data.push(format_cell(cell, headers[row_data.length]))
    })

    formatted_rows.push(row_data)
  })

  return { headers, rows: formatted_rows }
}

// Convert a table cell to its exported value: numbers (from data-sort-value or
// text) are formatted per the column's label spec; otherwise return cleaned text
function format_cell(cell: Element, header: string): string | number {
  // Prefer data-sort-value attribute (holds the unformatted number)
  const sort_value = cell.getAttribute(`data-sort-value`)
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

// Simple formatting function for export values
function format_number(value: number, format?: string): number | string {
  // SI/scientific formats use exponential notation for large numbers
  if ((format?.includes(`~s`) || format?.includes(`.3s`)) && Math.abs(value) >= 1000)
    return value.toExponential(2)

  // Parse decimal places from format strings like '.2f' (default 3)
  const n_decimals = Number(format?.match(/\.(?<decimals>\d)f/)?.groups?.decimals ?? 3)
  const factor = 10 ** n_decimals
  return Math.round(value * factor) / factor
}

// Helper function to format values according to label specifications
function format_value_for_export(value: number, header: string): number | string {
  // Find the corresponding label configuration by header text
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

    // Convert to CSV format
    const csv_data = [headers, ...rows]
    const csv_content = csv_data
      .map((row) =>
        row
          .map((cell) => {
            const cell_str = String(cell)
            // Escape quotes and wrap in quotes if necessary
            if (
              cell_str.includes(`,`) ||
              cell_str.includes(`"`) ||
              cell_str.includes(`\n`)
            ) {
              return `"${cell_str.replaceAll(`"`, `""`)}"`
            }
            return cell_str
          })
          .join(`,`),
      )
      .join(`\n`)

    // Create and download CSV file
    const blob = new Blob([csv_content], { type: `text/csv;charset=utf-8;` })
    const filename = generate_filename(`csv`, discovery_set)
    return download_blob(blob, filename)
  } catch (error) {
    console.error(`Error generating CSV:`, error)
    return null
  }
}

export async function generate_excel({
  discovery_set = `unique_prototypes`,
}: ExportOptions): Promise<ExportResult | null> {
  try {
    // Dynamic import of xlsx library
    const XLSX = await import(`xlsx`)

    const { headers, rows } = extract_table_data()

    // Create workbook and worksheet
    const excel_data = [headers, ...rows]
    const workbook = XLSX.utils.book_new()
    const worksheet = XLSX.utils.aoa_to_sheet(excel_data)

    // Auto-fit column widths (only if worksheet has data)
    if (worksheet[`!ref`]) {
      const range = XLSX.utils.decode_range(worksheet[`!ref`])
      const column_widths: { wch: number }[] = []

      for (let col = range.s.c; col <= range.e.c; col++) {
        let max_width = 10 // Minimum width

        for (let row = range.s.r; row <= range.e.r; row++) {
          const cell_address = XLSX.utils.encode_cell({ r: row, c: col })
          const cell = worksheet[cell_address]

          if (cell?.v != null) {
            const cell_length = String(cell.v).length
            max_width = Math.max(max_width, Math.min(cell_length, 50)) // Cap at 50 chars
          }
        }

        column_widths.push({ wch: max_width })
      }

      worksheet[`!cols`] = column_widths
    }

    // Add worksheet to workbook
    XLSX.utils.book_append_sheet(workbook, worksheet, `Metrics Table`)

    // Generate Excel file
    const excel_buffer = XLSX.write(workbook, { bookType: `xlsx`, type: `array` })
    const blob = new Blob([excel_buffer], {
      type: `application/vnd.openxmlformats-officedocument.spreadsheetml.sheet`,
    })
    const filename = generate_filename(`xlsx`, discovery_set)
    return download_blob(blob, filename)
  } catch (error) {
    console.error(`Error generating Excel:`, error)
    return null
  }
}

export const handle_export =
  <T extends ExportOptions>(
    generator: (args: T) => ExportResult | null | Promise<ExportResult | null>,
    fmt: string,
    state: { export_error: string | null } & T,
  ) =>
  async () => {
    try {
      state.export_error = null // Reset error state before trying
      // Extract only ExportOptions keys — T extends ExportOptions so this is safe
      const generator_args = {
        discovery_set: state.discovery_set,
      } as T // unavoidable: TS can't prove ExportOptions subset satisfies generic T
      const result = await generator(generator_args)
      if (!result) {
        state.export_error = `Failed to generate ${fmt}. The export function returned null.`
      }
    } catch (error) {
      state.export_error = `Error exporting ${fmt}: ${
        error instanceof Error ? error.message : String(error)
      }`
      console.error(`Error exporting ${fmt}:`, error)
    }
  }
