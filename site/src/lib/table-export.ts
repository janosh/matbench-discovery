import { toPng, toSvg } from 'html-to-image'
import { ALL_METRICS, HYPERPARAMS, METADATA_COLS } from './labels'

export const heatmap_class = `heatmap`

type ExportOptions = { show_non_compliant?: boolean; discovery_set?: string }

type ExportResult = { filename: string; url: string } | null

// Helper function to create a filtered table clone excluding SVG icon columns
function create_filtered_table_clone(): HTMLElement {
  const table_el = document.querySelector(`.${heatmap_class}`) as HTMLTableElement

  if (!table_el) throw new Error(`Table element not found for export`)

  // Clone the entire table
  const table_clone = table_el.cloneNode(true) as HTMLTableElement

  // Get header row to identify columns to remove
  const header_rows = table_clone.querySelectorAll(`thead tr`)
  const main_header_row = header_rows[header_rows.length - 1]
  const all_headers = Array.from(main_header_row.querySelectorAll(`th`))

  // Find column indices to remove (SVG icon columns)
  const columns_to_remove: number[] = []
  all_headers.forEach((th, index) => {
    const header_text = th.textContent?.replace(/[↑↓]/g, ``).trim() || ``
    if (header_text === `Org` || header_text === `Links`) {
      columns_to_remove.push(index)
    }
  })

  // Remove columns from all header rows (working backwards to maintain indices)
  columns_to_remove.reverse().forEach((col_index) => {
    header_rows.forEach((row) => {
      const cells = row.querySelectorAll(`th`)
      if (cells[col_index]) {
        cells[col_index].remove()
      }
    })
  })

  // Remove columns from all body rows
  const body_rows = table_clone.querySelectorAll(`tbody tr`)
  columns_to_remove.forEach((col_index) => {
    body_rows.forEach((row) => {
      const cells = row.querySelectorAll(`td`)
      if (cells[col_index]) {
        cells[col_index].remove()
      }
    })
  })

  return table_clone
}

// Remove HTML comments recursively (Svelte generates many of these)
function remove_comments(node: Node): void {
  if (node.nodeType === Node.COMMENT_NODE) {
    node.parentNode?.removeChild(node)
    return
  }

  // Recursively process child nodes (make a copy of the list since we're modifying it)
  const children = Array.from(node.childNodes)
  children.forEach(remove_comments)
}

// Clean up Svelte-specific artifacts and problematic elements
function clean_table_for_export(table_clone: HTMLElement): void {
  // Remove problematic elements
  const problematic_elements = table_clone.querySelectorAll(`svg, img, button, a[href]`)
  problematic_elements.forEach((el) => {
    // Replace with text content if available, otherwise remove entirely
    if (el.textContent?.trim()) {
      const text_node = document.createTextNode(el.textContent.trim())
      el.parentNode?.replaceChild(text_node, el)
    } else {
      el.remove()
    }
  })

  // Remove HTML comments
  remove_comments(table_clone)

  // Remove inline elements that might cause issues, but preserve sub/sup for visual formatting
  const inline_elements = table_clone.querySelectorAll(`span, small`)
  inline_elements.forEach((el) => {
    if (el.textContent?.trim()) {
      const text_node = document.createTextNode(el.textContent.trim())
      el.parentNode?.replaceChild(text_node, el)
    } else {
      el.remove()
    }
  })

  // Simplify class names and remove data attributes
  const all_elements = [table_clone, ...Array.from(table_clone.querySelectorAll(`*`))]
  all_elements.forEach((el) => {
    // Remove Svelte-generated dynamic classes
    if (el.className) {
      const classes = el.className.split(` `).filter((cls) => !cls.startsWith(`s-`))
      el.className = classes.join(` `)
    }

    // Remove Svelte-specific data attributes
    Array.from(el.attributes).forEach((attr) => {
      if (
        attr.name.startsWith(`data-`) &&
        ![`data-sort-value`, `data-title`].includes(attr.name)
      ) {
        el.removeAttribute(attr.name)
      }
    })

    // Clean up problematic styling
    const html_el = el as HTMLElement
    if (html_el.style) {
      html_el.style.removeProperty(`grid-column`)
      html_el.style.removeProperty(`grid-row`)
      html_el.style.removeProperty(`grid-area`)
    }
  })
}

// Create export container with proper styling
function create_export_container(table_clone: HTMLElement): HTMLElement {
  const night_color =
    getComputedStyle(document.documentElement).getPropertyValue(`--night`).trim() ||
    `#061e25`

  const container = document.createElement(`div`)
  container.style.backgroundColor = night_color
  container.style.display = `inline-block`
  container.style.whiteSpace = `nowrap`

  const inner = document.createElement(`div`)
  inner.style.padding = `15px`
  inner.style.boxSizing = `border-box`
  inner.appendChild(table_clone)
  container.appendChild(inner)

  return container
}

// Generate filename for export
function generate_filename(
  format: string,
  show_non_compliant: boolean = false,
  discovery_set: string = `unique_prototypes`,
): string {
  const date = new Date().toISOString().split(`T`)[0]
  const suffix = show_non_compliant ? `` : `-only-compliant`
  const param_case_discovery_set = discovery_set.replaceAll(`_`, `-`)
  return `metrics-table-${param_case_discovery_set}${suffix}-${date}.${format.toLowerCase()}`
}

// Common filter function for html-to-image
function create_export_filter() {
  return (node: Node) => {
    const node_name = node.nodeName

    // Skip external stylesheets
    if (
      node_name === `LINK` &&
      (node as Element).getAttribute(`rel`) === `stylesheet` &&
      (node as Element).getAttribute(`href`)?.startsWith(`http`)
    ) {
      return false
    }

    // Skip problematic elements
    if ([`SVG`, `IMG`].includes(node_name)) {
      return false
    }

    return true
  }
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

// Function to generate SVG from the metrics table
export async function generate_svg({
  show_non_compliant = false,
  discovery_set = `unique_prototypes`,
}: ExportOptions): Promise<ExportResult> {
  try {
    // Find the metrics table
    const table_el = document.querySelector(`.${heatmap_class}`)
    if (!table_el) {
      console.error(`Table element not found for SVG export`)
      return null
    }

    // Create and clean table clone
    const table_clone = create_filtered_table_clone()
    clean_table_for_export(table_clone)

    // Create export container
    const container = create_export_container(table_clone)
    document.body.appendChild(container)

    try {
      // Generate SVG
      const svg_data_url = await toSvg(container, {
        backgroundColor: container.style.backgroundColor,
        skipFonts: true,
        filter: create_export_filter(),
      })

      // Create download
      const blob = await fetch(svg_data_url).then((res) => res.blob())
      const url = URL.createObjectURL(blob)
      const filename = generate_filename(`svg`, show_non_compliant, discovery_set)

      const a = document.createElement(`a`)
      a.href = url
      a.download = filename
      a.click()

      setTimeout(() => URL.revokeObjectURL(url), 100)
      return { filename, url }
    } catch (error) {
      log_export_error(error, `SVG`, {
        containerHTML: container.outerHTML.substring(0, 1000) + `...`,
        tableCloneHTML: table_clone.outerHTML.substring(0, 1000) + `...`,
      })
      return null
    } finally {
      if (document.body.contains(container)) {
        document.body.removeChild(container)
      }
    }
  } catch (error) {
    log_export_error(error, `SVG setup`)
    return null
  }
}

// Function to generate PNG from the metrics table
export async function generate_png({
  show_non_compliant = false,
  discovery_set = `unique_prototypes`,
}: ExportOptions): Promise<ExportResult> {
  try {
    // Create and clean table clone
    const table_clone = create_filtered_table_clone()
    clean_table_for_export(table_clone)

    // Create export container with PNG-specific styling
    const container = create_export_container(table_clone)
    container.style.fontFamily = `Arial, sans-serif`
    container.style.fontSize = `14px`
    container.style.color = `#ffffff`

    document.body.appendChild(container)

    // Calculate dimensions for high-quality output
    const table_rect = table_clone.getBoundingClientRect()
    const padding = 15
    const width_with_padding = table_rect.width + padding * 2
    const height_with_padding = table_rect.height + padding * 2

    container.style.width = `${width_with_padding}px`
    container.style.height = `${height_with_padding}px`

    try {
      // Generate PNG
      const png_data_url = await toPng(container, {
        backgroundColor: container.style.backgroundColor,
        pixelRatio: 2, // High resolution
        width: width_with_padding,
        height: height_with_padding,
        skipFonts: true,
        quality: 0.95,
        filter: create_export_filter(),
        imagePlaceholder: `data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQDwAEhQGAhKmMIQAAAABJRU5ErkJggg==`,
      })

      // Create download
      const filename = generate_filename(`png`, show_non_compliant, discovery_set)
      const a = document.createElement(`a`)
      a.href = png_data_url
      a.download = filename
      a.click()

      return { filename, url: png_data_url }
    } catch (error) {
      log_export_error(error, `PNG`, {
        containerHTML: container.outerHTML.substring(0, 1000) + `...`,
        tableCloneHTML: table_clone.outerHTML.substring(0, 1000) + `...`,
      })
      return null
    } finally {
      if (document.body.contains(container)) {
        document.body.removeChild(container)
      }
    }
  } catch (error) {
    log_export_error(error, `PNG setup`)
    return null
  }
}

// Helper function to extract table data excluding SVG icon columns
function extract_table_data(): { headers: string[]; rows: (string | number)[][] } {
  const table_el = document.querySelector(`.${heatmap_class}`) as HTMLTableElement

  if (!table_el) {
    throw new Error(`Table element not found for export`)
  }

  // Extract headers from table
  const header_rows = table_el.querySelectorAll(`thead tr`)
  const main_header_row = header_rows[header_rows.length - 1]
  const all_headers = Array.from(main_header_row.querySelectorAll(`th`))

  // Filter out columns with SVG icons (Org and Links columns)
  const text_headers: string[] = []
  const included_column_indices: number[] = []

  all_headers.forEach((th, index) => {
    const header_text = th.textContent?.replace(/[↑↓]/g, ``).trim() || ``

    // Exclude columns that typically contain only SVG icons
    if (header_text !== `Org` && header_text !== `Links`) {
      text_headers.push(header_text)
      included_column_indices.push(index)
    }
  })

  // Extract data rows
  const data_rows = table_el.querySelectorAll(`tbody tr`)
  const formatted_rows: (string | number)[][] = []

  data_rows.forEach((row) => {
    const all_cells = Array.from(row.querySelectorAll(`td`))
    const row_data: (string | number)[] = []

    // Only include cells from columns we want to export
    included_column_indices.forEach((col_index) => {
      const cell = all_cells[col_index]
      const header = text_headers[row_data.length] // Get corresponding header

      if (cell) {
        let cell_value: string | number = ``

        // Check if cell has a data-sort-value attribute (for formatted numbers)
        const sort_value = cell.getAttribute(`data-sort-value`)
        if (sort_value && sort_value !== `null`) {
          // Try to parse as number
          const num_value = Number(sort_value)
          if (!isNaN(num_value)) {
            // Apply formatting based on column type
            cell_value = format_value_for_export(num_value, header)
          } else {
            cell_value = sort_value
          }
        } else {
          // Extract text content, handling HTML content
          const text_content = cell.textContent?.trim() || ``

          // Try to parse as number if it looks like one
          const num_value = Number(text_content)
          if (!isNaN(num_value) && text_content !== `` && text_content !== `n/a`) {
            cell_value = format_value_for_export(num_value, header)
          } else {
            cell_value = text_content
              .replace(/&lt;/g, `<`)
              .replace(/&gt;/g, `>`)
              .replace(/&amp;/g, `&`)
              .replace(/\s+/g, ` `)
              .trim()
          }
        }

        row_data.push(cell_value)
      }
    })

    formatted_rows.push(row_data)
  })

  return { headers: text_headers, rows: formatted_rows }
}

// Simple formatting function for export values
function format_number(value: number, format?: string): number | string {
  if (!format) {
    // Default: round to 3 decimal places
    return Math.round(value * 1000) / 1000
  }

  // Parse basic format strings
  if (format.includes(`.3f`)) {
    return Math.round(value * 1000) / 1000
  } else if (format.includes(`.2f`)) {
    return Math.round(value * 100) / 100
  } else if (format.includes(`.1f`)) {
    return Math.round(value * 10) / 10
  } else if (format.includes(`~s`) || format.includes(`.3s`)) {
    // Scientific notation for large numbers
    if (Math.abs(value) >= 1000) {
      return value.toExponential(2)
    }
    return Math.round(value * 1000) / 1000
  }

  // Default fallback
  return Math.round(value * 1000) / 1000
}

// Helper function to format values according to label specifications
function format_value_for_export(value: number, header: string): number | string {
  // Find the corresponding label configuration
  const all_labels = { ...ALL_METRICS, ...METADATA_COLS, ...HYPERPARAMS }

  // Look for label by header text (may need to handle HTML in headers)
  const clean_header = header.replace(/<[^>]*>/g, ``).trim()

  let format_spec: string | undefined

  // Find matching label by key, label, or short name
  for (const label of Object.values(all_labels)) {
    const label_text = label.label?.replace(/<[^>]*>/g, ``).trim()
    const short_text = label.short?.replace(/<[^>]*>/g, ``).trim()

    if (
      label_text === clean_header ||
      short_text === clean_header ||
      label.key === clean_header
    ) {
      format_spec = label.format
      break
    }
  }

  // Apply formatting
  return format_number(value, format_spec)
}

// Function to export table data as CSV
export async function generate_csv({
  show_non_compliant = false,
  discovery_set = `unique_prototypes`,
}: ExportOptions): Promise<ExportResult> {
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
              return `"${cell_str.replace(/"/g, `""`)}"`
            }
            return cell_str
          })
          .join(`,`),
      )
      .join(`\n`)

    // Create and download CSV file
    const blob = new Blob([csv_content], { type: `text/csv;charset=utf-8;` })
    const url = URL.createObjectURL(blob)
    const filename = generate_filename(`csv`, show_non_compliant, discovery_set)

    // Download the CSV
    const a = document.createElement(`a`)
    a.href = url
    a.download = filename
    a.click()

    // Clean up
    setTimeout(() => URL.revokeObjectURL(url), 100)

    return { filename, url }
  } catch (error) {
    console.error(`Error generating CSV:`, error)
    return null
  }
}

// Function to export table data as Excel (XLSX)
export async function generate_excel({
  show_non_compliant = false,
  discovery_set = `unique_prototypes`,
}: ExportOptions): Promise<ExportResult> {
  try {
    // Dynamic import of xlsx library
    const XLSX = await import(`xlsx`)

    const { headers, rows } = extract_table_data()

    // Create workbook and worksheet
    const excel_data = [headers, ...rows]
    const workbook = XLSX.utils.book_new()
    const worksheet = XLSX.utils.aoa_to_sheet(excel_data)

    // Auto-fit column widths
    const range = XLSX.utils.decode_range(worksheet[`!ref`] || `A1`)
    const column_widths: { wch: number }[] = []

    for (let col = range.s.c; col <= range.e.c; col++) {
      let max_width = 10 // minimum width

      for (let row = range.s.r; row <= range.e.r; row++) {
        const cell_address = XLSX.utils.encode_cell({ r: row, c: col })
        const cell = worksheet[cell_address]

        if (cell && cell.v) {
          const cell_length = String(cell.v).length
          max_width = Math.max(max_width, Math.min(cell_length, 50)) // cap at 50 chars
        }
      }

      column_widths.push({ wch: max_width })
    }

    worksheet[`!cols`] = column_widths

    // Add worksheet to workbook
    const sheet_name = `Metrics Table`
    XLSX.utils.book_append_sheet(workbook, worksheet, sheet_name)

    // Generate Excel file
    const excel_buffer = XLSX.write(workbook, { bookType: `xlsx`, type: `array` })
    const blob = new Blob([excel_buffer], {
      type: `application/vnd.openxmlformats-officedocument.spreadsheetml.sheet`,
    })
    const url = URL.createObjectURL(blob)
    const filename = generate_filename(`xlsx`, show_non_compliant, discovery_set)

    // Download the Excel file
    const a = document.createElement(`a`)
    a.href = url
    a.download = filename
    a.click()

    // Clean up
    setTimeout(() => URL.revokeObjectURL(url), 100)

    return { filename, url }
  } catch (error) {
    console.error(`Error generating Excel:`, error)
    return null
  }
}

export const handle_export =
  <T extends ExportOptions>(
    generator: (args: T) => Promise<ExportResult>,
    fmt: string,
    state: { export_error: string | null } & T,
  ) =>
  async () => {
    try {
      state.export_error = null // Reset error state before trying
      // Pass only the relevant properties expected by the generator
      const generator_args: T = {
        show_non_compliant: state.show_non_compliant,
        discovery_set: state.discovery_set,
      } as T // Cast needed as state has extra key
      const result = await generator(generator_args)
      if (!result) {
        state.export_error = `Failed to generate ${fmt}. The export function returned null.`
      }
    } catch (err) {
      state.export_error = `Error exporting ${fmt}: ${err instanceof Error ? err.message : String(err)}`
      console.error(`Error exporting ${fmt}:`, err)
    }
  }
