import { toPng, toSvg } from 'html-to-image'

// Function to export as SVG using html-to-image
export async function generate_svg({
  show_non_compliant = false,
  discovery_set = `unique_prototypes`,
}: {
  show_non_compliant?: boolean
  discovery_set?: string
}) {
  try {
    // Find the metrics table
    const table_el = document.querySelector(`.heatmap-table`)

    if (!table_el) {
      console.error(`Table element not found for SVG export`)
      return null
    }

    // Get the night theme color for the background
    const night_color = getComputedStyle(document.documentElement)
      .getPropertyValue(`--night`) // see app.css
      .trim()

    // Clone the table first so we can measure it properly
    const table_clone = table_el.cloneNode(true) as HTMLElement

    // Create a simple container with the table clone
    const container = document.createElement(`div`)
    container.style.backgroundColor = night_color
    container.style.display = `inline-block` // This makes container size to content
    container.style.whiteSpace = `nowrap` // Prevent text wrapping

    // Create an inner container with padding to ensure even margins
    const inner = document.createElement(`div`)
    inner.style.padding = `15px` // Consistent padding on all sides
    inner.style.boxSizing = `border-box` // Include padding in width/height
    inner.appendChild(table_clone)
    container.appendChild(inner)

    // Temporarily add to document to calculate dimensions
    document.body.appendChild(container)

    // Get the table dimensions
    const table_rect = table_clone.getBoundingClientRect()

    // Add padding to table dimensions - 15px on each side
    const padding = 15 // Must match padding in inner div
    const table_width_with_padding = table_rect.width + padding * 2
    const table_height_with_padding = table_rect.height + padding * 2

    // Set container dimensions based on table size plus padding
    container.style.width = `${table_width_with_padding}px`
    container.style.height = `${table_height_with_padding}px`

    try {
      // Generate SVG with html-to-image
      const svg_data_url = await toSvg(container, {
        backgroundColor: night_color,
        width: table_width_with_padding,
        height: table_height_with_padding,
        skipFonts: true, // Skip remote font processing
        filter: (node) => {
          // Skip processing of link elements with external stylesheets
          if (
            node.nodeName === `LINK` &&
            node.getAttribute(`rel`) === `stylesheet` &&
            node.getAttribute(`href`)?.startsWith(`http`)
          ) {
            return false
          }
          return true
        },
      })

      // SVG data URL to Blob - using fetch to avoid encoding issues
      const blob = await fetch(svg_data_url).then((res) => res.blob())
      const url = URL.createObjectURL(blob)

      // Create filename with date
      const date = new Date().toISOString().split(`T`)[0]
      const suffix = show_non_compliant ? `` : `-only-compliant`
      const param_case_discovery_set = discovery_set.replaceAll(`_`, `-`)
      const filename = `metrics-table-${param_case_discovery_set}${suffix}-${date}.svg`

      // Download the SVG
      const a = document.createElement(`a`)
      a.href = url
      a.download = filename
      a.click()

      // Clean up
      setTimeout(() => URL.revokeObjectURL(url), 100)

      return { filename, url }
    } catch (error) {
      console.error(`Error exporting SVG:`, error)
      return null
    } finally {
      // Always clean up the container
      if (document.body.contains(container)) {
        document.body.removeChild(container)
      }
    }
  } catch (error) {
    console.error(`Error generating SVG:`, error)
    return null
  }
}

// Function to generate a high-resolution PNG from the metrics table
export async function generate_png({
  show_non_compliant = false,
  discovery_set = `unique_prototypes`,
}: {
  show_non_compliant?: boolean
  discovery_set?: string
}) {
  try {
    // Find the metrics table
    const table_el = document.querySelector(`.heatmap-table`)

    if (!table_el) {
      console.error(`Table element not found for PNG export`)
      return null
    }

    // Get the night theme color for the background
    const night_color =
      getComputedStyle(document.documentElement).getPropertyValue(`--night`).trim() ||
      `#061e25`

    // Clone the table first so we can measure it properly
    const table_clone = table_el.cloneNode(true) as HTMLElement

    // Create a simple container with the table clone
    const container = document.createElement(`div`)
    container.style.backgroundColor = night_color
    container.style.display = `inline-block` // This makes container size to content
    container.style.whiteSpace = `nowrap` // Prevent text wrapping

    // Create an inner container with padding to ensure even margins
    const inner = document.createElement(`div`)
    inner.style.padding = `15px` // Consistent padding on all sides
    inner.style.boxSizing = `border-box` // Include padding in width/height
    inner.appendChild(table_clone)
    container.appendChild(inner)

    // Temporarily add to document to calculate dimensions
    document.body.appendChild(container)

    // Get the table dimensions first
    const table_rect = table_clone.getBoundingClientRect()

    // Add padding to table dimensions - 15px on each side
    const padding = 15 // Must match padding in inner div
    const table_width_with_padding = table_rect.width + padding * 2
    const table_height_with_padding = table_rect.height + padding * 2

    // Set container dimensions based on table size plus padding
    container.style.width = `${table_width_with_padding}px`
    container.style.height = `${table_height_with_padding}px`

    try {
      // Generate PNG with html-to-image
      const png_data_url = await toPng(container, {
        backgroundColor: night_color,
        pixelRatio: 2, // High resolution (2x)
        width: table_width_with_padding,
        height: table_height_with_padding,
        skipFonts: true, // Skip remote font processing
        filter: (node) => {
          // Skip processing of link elements with external stylesheets
          if (
            node.nodeName === `LINK` &&
            node.getAttribute(`rel`) === `stylesheet` &&
            node.getAttribute(`href`)?.startsWith(`http`)
          ) {
            return false
          }
          return true
        },
        imagePlaceholder: `data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQDwAEhQGAhKmMIQAAAABJRU5ErkJggg==`,
      })

      // Create filename with date
      const date = new Date().toISOString().split(`T`)[0]
      const suffix = show_non_compliant ? `` : `-only-compliant`
      const param_case_discovery_set = discovery_set.replaceAll(`_`, `-`)
      const filename = `metrics-table-${param_case_discovery_set}${suffix}-${date}.png`

      // Download the PNG
      const a = document.createElement(`a`)
      a.href = png_data_url
      a.download = filename
      a.click()

      return { filename, url: png_data_url }
    } catch (error) {
      console.error(`Error exporting PNG:`, error)
      return null
    } finally {
      // Always clean up the container
      if (document.body.contains(container)) {
        document.body.removeChild(container)
      }
    }
  } catch (error) {
    console.error(`Error generating PNG:`, error)
    return null
  }
}

export const handle_export =
  (
    generator: typeof generate_svg | typeof generate_png,
    fmt: string,
    export_error: string | null,
    state: Record<string, unknown>,
  ) =>
  async () => {
    try {
      export_error = null
      const result = await generator(state)
      if (!result) {
        export_error = `Failed to generate ${fmt}. The export function returned null.`
      }
    } catch (error) {
      error = `Error exporting ${fmt}: ${error instanceof Error ? error.message : String(error)}`
      console.error(`Error exporting ${fmt}:`, error)
    }
  }
