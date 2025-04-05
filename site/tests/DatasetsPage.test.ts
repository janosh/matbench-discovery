import { mount, tick } from 'svelte'
import { beforeEach, describe, expect, it } from 'vitest'
import { doc_query } from '.'
import Page from '../src/routes/data/sets/+page.svelte'

describe(`Datasets Page`, () => {
  beforeEach(() => {
    mount(Page, { target: document.body })
  })

  it(`renders the page title correctly`, () => {
    const heading = doc_query<HTMLHeadingElement>(`h1`)
    expect(heading.textContent).toBe(`Datasets`)
  })

  it(`renders the table with correct structure`, () => {
    const table = doc_query<HTMLTableElement>(`.training-sets-table`)
    const thead = table.querySelector(`thead`)
    const tbody = table.querySelector(`tbody`)

    // Check table structure
    expect(thead).not.toBeNull()
    expect(tbody).not.toBeNull()

    // Check header columns
    const header_cols = thead?.querySelectorAll(`th`) || []
    expect(header_cols.length).toBe(9)

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
    const rows = document.querySelectorAll(`tbody tr`)
    expect(rows.length).toBeGreaterThan(5)

    // Check for some expected dataset names in the table
    const dataset_names = Array.from(rows).map((row) => {
      const title_cell = row.querySelector(`td[data-col="Title"]`)
      return title_cell?.textContent?.trim() || ``
    })

    expect(dataset_names).toContain(`MP 2022`)
    expect(dataset_names.some((name) => name.includes(`Alex`))).toBe(true)
  })

  it(`properly renders links for datasets`, () => {
    // Check that all dataset rows have link cells
    const link_cells = document.querySelectorAll(`td[data-col="Links"]`)
    expect(link_cells.length).toBeGreaterThan(5)

    // Verify links are present and have correct structure
    let link_count = 0
    link_cells.forEach((cell) => {
      const links = cell.querySelectorAll(`a`)
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
    // Get initial order of datasets by title
    const initial_datasets = Array.from(document.querySelectorAll(`tbody tr`)).map(
      (row) => {
        const title_cell = row.querySelector(`td[data-col="Title"]`)
        return title_cell?.textContent?.trim() || ``
      },
    )

    // Find and click the Structures column header to sort
    const structures_header = Array.from(document.querySelectorAll(`th`)).find((th) =>
      th.textContent?.includes(`Structures`),
    )
    structures_header?.click()
    await tick()

    // Get new order after sorting
    const sorted_datasets = Array.from(document.querySelectorAll(`tbody tr`)).map(
      (row) => {
        const title_cell = row.querySelector(`td[data-col="Title"]`)
        return title_cell?.textContent?.trim() || ``
      },
    )

    // Order should have changed
    expect(sorted_datasets).not.toEqual(initial_datasets)
  })

  it(`skips sorting when clicking non-sortable columns`, async () => {
    // Get initial order of datasets
    const initial_datasets = Array.from(document.querySelectorAll(`tbody tr`)).map(
      (row) => {
        const title_cell = row.querySelector(`td[data-col="Title"]`)
        return title_cell?.textContent?.trim() || ``
      },
    )

    // Find and click the Links column header (which should be non-sortable)
    const links_header = Array.from(document.querySelectorAll(`th`)).find((th) =>
      th.textContent?.includes(`Links`),
    )
    links_header?.click()
    await tick()

    // Get new order after clicking non-sortable column
    const new_datasets = Array.from(document.querySelectorAll(`tbody tr`)).map((row) => {
      const title_cell = row.querySelector(`td[data-col="Title"]`)
      return title_cell?.textContent?.trim() || ``
    })

    // Order should remain the same
    expect(new_datasets).toEqual(initial_datasets)
  })

  it(`has correct styling for sortable and non-sortable columns`, () => {
    const sortable_headers = document.querySelectorAll(`th.sortable`)
    const all_headers = document.querySelectorAll(`th`)

    // Most columns should be sortable
    expect(sortable_headers.length).toBe(all_headers.length - 1)

    // Links column should not have sortable class
    const links_header = Array.from(all_headers).find((th) =>
      th.textContent?.includes(`Links`),
    )
    expect(links_header?.classList.contains(`sortable`)).toBe(false)
  })

  it(`formats numbers correctly in the table`, () => {
    const structures_cells = document.querySelectorAll(`td[data-col="Structures"]`)

    // Check that at least some cells have formatted numbers
    let has_formatted_number = false
    structures_cells.forEach((cell) => {
      const cell_text = cell.textContent?.trim() || ``
      if (cell_text !== `n/a`) {
        // Should use K or M for thousands/millions
        has_formatted_number = has_formatted_number || /\d+(\.\d+)?[KM]/.test(cell_text)
      }
    })

    expect(has_formatted_number).toBe(true)
  })
})
