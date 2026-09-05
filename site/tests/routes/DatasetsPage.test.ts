import { DATASETS } from '$lib'
import { heatmap_class } from '$lib/table-export'
import Page from '$routes/data/sets/+page.svelte'
import { beforeEach, describe, expect, it } from 'vitest'
import { doc_query, mount } from '../index'

describe(`Datasets Page`, () => {
  beforeEach(() => {
    mount(Page, { target: document.body })
  })

  it(`renders the table with correct structure`, () => {
    expect(doc_query<HTMLHeadingElement>(`h1`).textContent).toContain(`Datasets`)
    const table = doc_query<HTMLTableElement>(`.${heatmap_class}`)
    const thead = doc_query(`thead`, table)
    const tbody = doc_query(`tbody`, table)

    const header_cols = thead.querySelectorAll(`th`)
    expect(header_cols).toHaveLength(10)

    const column_headers = [...header_cols].map(
      (col) => col.textContent?.trim().replaceAll(/[↑↓]/g, ``) ?? ``,
    )
    expect(column_headers.some((header) => header.includes(`Name`))).toBe(true)
    expect(column_headers.some((header) => header.includes(`Structures`))).toBe(true)
    expect(column_headers.some((header) => header.includes(`Materials`))).toBe(true)
    expect(column_headers.some((header) => header.includes(`Created`))).toBe(true)
    expect(column_headers.some((header) => header.includes(`API`))).toBe(true)
    expect(column_headers.some((header) => header.includes(`Links`))).toBe(true)

    // One row per dataset in datasets.yml, title in the first cell
    const rows = tbody.querySelectorAll(`tr`)
    expect(rows).toHaveLength(Object.keys(DATASETS).length)
    const dataset_names = [...rows].map(
      (row) => row.querySelector(`td:first-child`)?.textContent?.trim() ?? ``,
    )
    expect(dataset_names.some((name) => name.includes(`MP`))).toBe(true)
    expect(dataset_names.some((name) => name.includes(`Alex`))).toBe(true)
  })

  it(`properly renders resource links for datasets`, () => {
    // Links column is the last column (10th, index 9)
    const tbody = document.querySelector(`.${heatmap_class} tbody`)
    if (!tbody) throw new Error(`Datasets table body not found`)
    const rows = tbody.querySelectorAll(`tr`)

    // Count resource links (Website, Download, DOI)
    let resource_link_count = 0
    rows.forEach((row) => {
      const cells = row.querySelectorAll(`td`)
      const links_cell = cells[9]
      if (!links_cell) return

      const resource_links = [...links_cell.querySelectorAll(`a`)].filter(
        (link) => link.hasAttribute(`title`) || link.hasAttribute(`data-original-title`),
      )
      resource_link_count += resource_links.length

      resource_links.forEach((link) => {
        expect(link.getAttribute(`target`)).toBe(`_blank`)
        expect(link.getAttribute(`rel`)).toContain(`noopener`)
        const title =
          link.getAttribute(`title`) ?? link.getAttribute(`data-original-title`)
        expect([`Website`, `Download`, `DOI`]).toContain(title)
        expect(link.getAttribute(`aria-label`)).toBe(title)
      })
    })

    expect(resource_link_count).toBeGreaterThan(10)
  })

  it(`properly renders API links for datasets`, () => {
    // API column is the 9th column (index 8)
    const rows = document.querySelectorAll(`.${heatmap_class} tbody tr`)
    const api_links = [...rows].flatMap((row) =>
      [...row.querySelectorAll(`td`)[8].children].filter(
        (child): child is HTMLAnchorElement => child.tagName === `A`,
      ),
    )

    // One API link per dataset native_api/optimade_api URL
    const by_string = (str_1: string, str_2: string) => str_1.localeCompare(str_2)
    const expected_api_hrefs = Object.values(DATASETS)
      .flatMap((dataset) => [dataset.native_api, dataset.optimade_api])
      .filter((href): href is string => Boolean(href))
      .toSorted(by_string)
    // raw attribute, since `.href` normalizes bare origins with a trailing slash
    const api_hrefs = api_links.map((link) => link.getAttribute(`href`) ?? ``)
    expect(api_hrefs.toSorted(by_string)).toStrictEqual(expected_api_hrefs)

    for (const link of api_links) {
      expect(link.getAttribute(`target`)).toBe(`_blank`)
      expect(link.getAttribute(`rel`)).toBe(`noopener noreferrer`)
      expect([`Native API`, `OPTIMADE API`]).toContain(link.getAttribute(`title`))
    }
  })

  it(`marks only the Links and API columns as non-sortable`, () => {
    // In HeatmapTable, non-sortable columns have the 'not-sortable' class
    const all_headers = [...document.querySelectorAll(`.${heatmap_class} th`)]
    const non_sortable = all_headers.filter((th) => th.classList.contains(`not-sortable`))
    expect(non_sortable.map((th) => th.textContent?.trim())).toStrictEqual([
      `API`,
      `Links`,
    ])
  })

  it(`formats numbers correctly in the table`, () => {
    // Find Structures column (usually 2nd column) cells
    const tbody = document.querySelector(`.${heatmap_class} tbody`)
    const rows = tbody?.querySelectorAll(`tr`) ?? []

    const structures_cells = [...rows].map((row) => {
      const cells = row.querySelectorAll(`td`)
      return cells[1]
    })

    // K/M suffixes for thousands/millions
    const has_formatted_number = structures_cells.some((cell) => {
      const cell_text = cell?.textContent?.trim() ?? ``
      return cell_text !== `n/a` && /\d+(?:\.\d+)?[KM]/.test(cell_text)
    })

    expect(has_formatted_number).toBe(true)
  })

  it(`correctly displays method information in the table`, () => {
    // Method is usually the 8th column (index 7)
    const tbody = document.querySelector(`.${heatmap_class} tbody`)
    const rows = tbody?.querySelectorAll(`tr`) ?? []

    const method_cells = [...rows].map((row) => {
      const cells = row.querySelectorAll(`td`)
      return cells[7]
    })

    const method_count = method_cells.filter(
      (cell) => cell?.textContent?.trim() !== `n/a`,
    ).length

    expect(method_count).toBeGreaterThan(3)

    const all_methods_text = method_cells
      .map((cell) => cell?.textContent?.trim())
      .join(` `)

    expect(/DFT|ML|experiment|GW|DMFT|MD/.test(all_methods_text)).toBe(true)
  })
})
