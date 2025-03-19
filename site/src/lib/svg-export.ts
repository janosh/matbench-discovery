// Function to generate SVG from HTML table
export async function generate_svg_from_table(): Promise<string | null> {
  const table_el =
    document.querySelector(`.heatmap-table`) ||
    document.querySelector(`.metrics-table table`) ||
    document.querySelector(`.metrics-table .heatmap-table-container table`)

  if (!table_el) {
    console.error(`Table element not found`)
    return null
  }

  try {
    // Log the number of rows to make sure we're getting all of them
    const rows = table_el.querySelectorAll(`tr`)

    // Calculate dimensions by summing the height of each row
    let [total_height, max_width] = [0, 0]

    // Add some extra spacing to account for headers
    total_height += 30 // Add space for headers

    for (const row of rows) {
      const row_height = row.getBoundingClientRect().height
      total_height += row_height

      const row_width = row.getBoundingClientRect().width
      if (row_width > max_width) max_width = row_width
    }

    // Add extra SVG size to ensure we capture all table HTML
    total_height += 50
    max_width += 50

    // Clone the table to avoid modifying the displayed one
    const table_clone = table_el.cloneNode(true) as HTMLElement

    // Apply styling for better appearance and compact rows
    table_clone.style.fontFamily = `Arial, sans-serif`
    table_clone.style.borderCollapse = `collapse` // Remove gaps between cells
    table_clone.style.width = `${max_width}px`
    table_clone.style.maxWidth = `${max_width}px`

    // Create a style element for more specific styling
    const style = document.createElement(`style`)
    style.textContent = `
        table { border-collapse: collapse; width: 100%; }
        th, td {
          padding: 2px 4px !important;
          line-height: 1 !important;
          font-size: 90% !important;
          border: 1px solid #ccc;
        }
        tr { height: auto !important; }
      `

    // Process all cells to capture their background colors and add borders
    const cells = table_clone.querySelectorAll(`td, th`)
    cells.forEach((cell) => {
      const cell_element = cell as HTMLElement
      const computed_style = window.getComputedStyle(cell)

      // Make rows more compact with minimal padding and line height
      cell_element.style.padding = `1px 2px`
      cell_element.style.lineHeight = `1`
      cell_element.style.height = `auto`

      // Copy background color from computed style explicitly
      const bg_color = computed_style.backgroundColor
      if (bg_color && bg_color !== `rgba(0, 0, 0, 0)` && bg_color !== `transparent`) {
        cell_element.style.backgroundColor = bg_color
      }

      // Add thin borders to all cells
      cell_element.style.border = `1px solid #ccc`
    })

    // Create a container to hold the table with the style
    const container = document.createElement(`div`)
    container.appendChild(style)
    container.appendChild(table_clone)

    const padding = 10 // white space around the table
    // Create SVG with foreignObject to embed the table as HTML and embedded style
    return `
      <svg xmlns="http://www.w3.org/2000/svg" width="${max_width + padding * 2}" height="${total_height + padding * 2}" viewBox="0 0 ${max_width + padding * 2} ${total_height + padding * 2}">
        <rect width="100%" height="100%" fill="white"/>
        <foreignObject x="${padding}" y="${padding}" width="${max_width}" height="${total_height}">
          ${new XMLSerializer().serializeToString(container)}
        </foreignObject>
      </svg>
      `
  } catch (error) {
    console.error(`Error generating SVG:`, error)
    return null
  }
}

// Function to export as SVG (reusing our SVG generator)
export async function generate_svg({
  show_non_compliant = false,
  discovery_set = `unique_prototypes`,
}: {
  show_non_compliant?: boolean
  discovery_set?: string
}) {
  try {
    const svg_data = await generate_svg_from_table()
    if (!svg_data) {
      console.error(`Failed to generate SVG: table element not found`)
      return null
    }

    // Create blob and download
    const blob = new Blob([svg_data], { type: `image/svg+xml` })
    const url = URL.createObjectURL(blob)

    // Create filename with date
    const date = new Date().toISOString().split(`T`)[0] // YYYY-MM-DD format
    const suffix = show_non_compliant ? `` : `-only-compliant`
    // Convert snake_case to param-case for the filename
    const param_case_discovery_set = discovery_set.replaceAll(`_`, `-`)
    const filename = `metrics-table-${param_case_discovery_set}${suffix}-${date}.svg`

    const a = document.createElement(`a`)
    a.href = url
    a.download = filename
    a.click()

    // Clean up
    setTimeout(() => URL.revokeObjectURL(url), 100)

    return { filename, url }
  } catch (error) {
    console.error(`Error generating SVG:`, error)
    return null
  }
}

export async function copy_pdf_conversion_cmd() {
  try {
    // First trigger SVG download and get the filename
    const svg = await generate_svg({
      show_non_compliant: false,
      discovery_set: `unique_prototypes`,
    })

    if (!svg) {
      console.error(`Failed to generate SVG`)
      return null
    }

    const pdf_filename = svg.filename.replace(`.svg`, `.pdf`)
    const base_filename = svg.filename.replace(`.svg`, ``)

    // Create the command with brace expansion syntax and installation instructions
    const command = `# Install pdf2svg:
# Linux: sudo apt-get install pdf2svg
# macOS: brew install pdf2svg
pdf2svg ${base_filename}.{svg,pdf}`

    // Copy to clipboard
    await navigator.clipboard.writeText(command)

    // Show a temporary success message
    const button = document.getElementById(`pdf-command-btn`)
    if (button) {
      const original_text = button.textContent
      button.textContent = `✓ Copied!`
      setTimeout(() => {
        button.textContent = original_text
      }, 2000)
    }

    return { command, svg_filename: svg.filename, pdf_filename }
  } catch (error) {
    console.error(`Error copying command:`, error)
    return null
  }
}
