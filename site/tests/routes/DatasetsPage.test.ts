import Page from '$routes/data/sets/+page.svelte'
import { mount, tick } from 'svelte'
import { beforeEach, describe, expect, it } from 'vitest'
import { doc_query } from '..'

describe(`Datasets Page`, () => {
  beforeEach(() => {
    mount(Page, { target: document.body })
  })

  it(`renders the page title correctly`, () => {
    const heading = doc_query<HTMLHeadingElement>(`h1`)
    expect(heading.textContent).toBe(`Datasets`)
  })

  it(`renders the table with correct structure`, () => {
    const table = doc_query<HTMLTableElement>(`.heatmap`)
    const thead = table.querySelector(`thead`)
    const tbody = table.querySelector(`tbody`)

    // Check table structure
    expect(thead).not.toBeNull()
    expect(tbody).not.toBeNull()

    // Check header columns
    const header_cols = thead?.querySelectorAll(`th`) || []
    expect(header_cols.length).toBe(8)

    // Verify expected column headers are present
    const column_headers = Array.from(header_cols).map(
      (col) => col.textContent?.trim().replace(/[↑↓]/g, ``) || ``,
    )
    expect(column_headers.some((header) => header.includes(`Title`))).toBe(true)
    expect(column_headers.some((header) => header.includes(`Structures`))).toBe(true)
    expect(column_headers.some((header) => header.includes(`Materials`))).toBe(true)
    expect(column_headers.some((header) => header.includes(`Created`))).toBe(true)
    expect(column_headers.some((header) => header.includes(`Links`))).toBe(true)

    // Check that we have rows in the table
    const rows = tbody?.querySelectorAll(`tr`) || []
    expect(rows.length).toBeGreaterThan(0)
  })

  it(`displays information about multiple datasets`, () => {
    const rows = document.querySelectorAll(`.heatmap tbody tr`)
    expect(rows.length).toBeGreaterThan(5)

    // Check for some expected dataset names in the table
    const dataset_names = Array.from(rows).map((row) => {
      // Title is always the first cell
      const title_cell = row.querySelector(`td:first-child`)
      return title_cell?.textContent?.trim() || ``
    })

    expect(dataset_names.some((name) => name.includes(`MP`))).toBe(true)
    expect(dataset_names.some((name) => name.includes(`Alex`))).toBe(true)
  })

  it(`properly renders links for datasets`, () => {
    // Links column is the last column (8th)
    const tbody = document.querySelector(`.heatmap tbody`)
    const rows = tbody?.querySelectorAll(`tr`) || []

    // Count links by checking each row's last cell (Links column)
    let link_count = 0
    rows.forEach((row) => {
      const cells = row.querySelectorAll(`td`)
      const links_cell = cells[cells.length - 1]
      const links = links_cell.querySelectorAll(`a`)
      link_count += links.length

      // Each link should have target="_blank" for external opening
      links.forEach((link) => {
        expect(link.getAttribute(`target`)).toBe(`_blank`)
        expect(link.getAttribute(`rel`)).toContain(`noopener`)
      })
    })

    // There should be multiple links across all cells
    expect(link_count).toBeGreaterThan(10)
  })

  it(`can sort table by clicking column headers`, async () => {
    // Get initial order of datasets
    const initial_datasets = Array.from(
      document.querySelectorAll(`.heatmap tbody tr`),
    ).map((row) => {
      // First cell is Title column
      const cells = row.querySelectorAll(`td`)
      return cells[0]?.textContent?.trim() || ``
    })

    // Find and click the Structures column header to sort (usually the 2nd header)
    const headers = document.querySelectorAll(`.heatmap th`)
    const structures_header = Array.from(headers).find((th) =>
      th.textContent?.includes(`Structures`),
    )
    if (structures_header) {
      ;(structures_header as HTMLElement).click()
    }
    await tick()

    // Get new order after sorting
    const sorted_datasets = Array.from(
      document.querySelectorAll(`.heatmap tbody tr`),
    ).map((row) => {
      const cells = row.querySelectorAll(`td`)
      return cells[0]?.textContent?.trim() || ``
    })

    // Order should have changed
    expect(sorted_datasets).not.toEqual(initial_datasets)
  })

  it(`skips sorting when clicking non-sortable columns`, async () => {
    // Get initial order of datasets
    const initial_datasets = Array.from(
      document.querySelectorAll(`.heatmap tbody tr`),
    ).map((row) => {
      const cells = row.querySelectorAll(`td`)
      return cells[0]?.textContent?.trim() || ``
    })

    // Find and click the Links column header (which should be non-sortable)
    const links_header = Array.from(document.querySelectorAll(`.heatmap th`)).find((th) =>
      th.textContent?.includes(`Links`),
    )
    if (links_header) {
      ;(links_header as HTMLElement).click()
    }
    await tick()

    // Get new order after clicking non-sortable column
    const new_datasets = Array.from(document.querySelectorAll(`.heatmap tbody tr`)).map(
      (row) => {
        const cells = row.querySelectorAll(`td`)
        return cells[0]?.textContent?.trim() || ``
      },
    )

    // Order should remain the same
    expect(new_datasets).toEqual(initial_datasets)
  })

  it(`has correct styling for sortable and non-sortable columns`, () => {
    // In HeatmapTable, non-sortable columns have the 'not-sortable' class
    const non_sortable_headers = document.querySelectorAll(`.heatmap th.not-sortable`)
    const all_headers = document.querySelectorAll(`.heatmap th`)

    // Only Links column should have not-sortable class
    expect(non_sortable_headers.length).toBe(1)

    // The Links header should have the not-sortable class
    const links_header = Array.from(all_headers).find((th) =>
      th.textContent?.includes(`Links`),
    )
    expect(links_header?.classList.contains(`not-sortable`)).toBe(true)
  })

  it(`formats numbers correctly in the table`, () => {
    // Find Structures column (usually 2nd column) cells
    const tbody = document.querySelector(`.heatmap tbody`)
    const rows = tbody?.querySelectorAll(`tr`) || []

    // Get cells from second column (Structures)
    const structures_cells = Array.from(rows).map((row) => {
      const cells = row.querySelectorAll(`td`)
      return cells[1] // Second column (index 1)
    })

    // Check that at least some cells have formatted numbers
    let has_formatted_number = false
    structures_cells.forEach((cell) => {
      const cell_text = cell?.textContent?.trim() || ``
      if (cell_text !== `n/a`) {
        // Should use K or M for thousands/millions
        has_formatted_number = has_formatted_number || /\d+(\.\d+)?[KM]/.test(cell_text)
      }
    })

    expect(has_formatted_number).toBe(true)
  })

  it(`correctly displays method information in the table`, () => {
    // Method is usually the 7th column (index 6)
    const tbody = document.querySelector(`.heatmap tbody`)
    const rows = tbody?.querySelectorAll(`tr`) || []

    // Get cells from the Method column
    const method_cells = Array.from(rows).map((row) => {
      const cells = row.querySelectorAll(`td`)
      return cells[6] // Method column
    })

    // At least some cells should have method information (not all n/a)
    const method_count = method_cells.filter(
      (cell) => cell?.textContent?.trim() !== `n/a`,
    ).length

    // There should be at least several datasets with method information
    expect(method_count).toBeGreaterThan(3)

    // Check for common methods like DFT or ML in the column
    const all_methods_text = method_cells
      .map((cell) => cell?.textContent?.trim())
      .join(` `)

    // Should find at least one of these common methods
    expect(/DFT|ML|experiment/.test(all_methods_text)).toBe(true)
  })
})
