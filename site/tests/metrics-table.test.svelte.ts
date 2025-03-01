import MetricsTable from '$lib/MetricsTable.svelte'
import type { HeatmapColumn, ModelData } from '$lib/types'
import { mount, tick } from 'svelte'
import { describe, expect, it } from 'vitest'

describe(`MetricsTable`, () => {
  it(`renders with default props`, () => {
    mount(MetricsTable, {
      target: document.body,
      props: {
        discovery_set: `unique_prototypes`,
        col_filter: () => true,
      },
    })

    // Check table structure
    const table = document.body.querySelector(`table`)
    expect(table).toBeDefined()
    expect(table?.querySelector(`thead`)).toBeDefined()
    expect(table?.querySelector(`tbody`)).toBeDefined()

    // Check number of columns
    const headers = document.body.querySelectorAll(`th`)
    expect(headers.length).toBeGreaterThan(15) // Should have all metric and metadata columns

    // Check essential columns are present (ignoring sort indicators)
    const header_texts = [...headers].map((h) => h.textContent?.trim())
    const required_cols = [
      `Model`,
      `F1 ↑`,
      `DAF ↑`,
      `Training Set`,
      `Params`,
      `Targets`,
      `Links`,
    ]
    expect(header_texts).toEqual(expect.arrayContaining(required_cols))

    // Check number of rows
    const rows = document.body.querySelectorAll(`tbody tr`)
    expect(rows.length).toBeGreaterThan(7) // Should have at least one model
  })

  it(`toggles metadata columns`, async () => {
    const metadata_cols = [`Training Set`, `Params`, `Targets`, `Date Added`, `Links`]
    let col_filter = $state((_col: HeatmapColumn) => true) // show all columns initially
    mount(MetricsTable, {
      target: document.body,
      props: { discovery_set: `unique_prototypes`, col_filter },
    })

    // Check metadata columns are visible
    let header_texts = [...document.body.querySelectorAll(`th`)].map((h) =>
      h.textContent?.trim(),
    )
    expect(header_texts).toEqual(expect.arrayContaining(metadata_cols))

    // Hide metadata columns
    col_filter = (col) => !metadata_cols.includes(col.label)
    await tick()

    // Check metadata columns are hidden
    header_texts = [...document.body.querySelectorAll(`th`)].map((h) =>
      h.textContent?.trim(),
    )
    for (const col of metadata_cols) {
      expect(header_texts).not.toContain(col)
    }

    // Check metric columns are still visible
    const metric_cols = [`F1 ↑`, `DAF ↑`, `Prec ↑`, `Acc ↑`]
    for (const col of metric_cols) {
      expect(header_texts).toContain(col)
    }
  })

  it(`filters specified columns`, async () => {
    mount(MetricsTable, {
      target: document.body,
      props: {
        discovery_set: `unique_prototypes`,
        col_filter: (col) => ![`F1`, `DAF`].includes(col.label),
      },
    })

    const headers = document.body.querySelectorAll(`th`)
    const header_texts = Array.from(headers).map((h) => h.textContent?.split(` `)[0])

    // Check hidden columns
    expect(header_texts).not.toContain(`F1`)
    expect(header_texts).not.toContain(`DAF`)

    // Check other columns still visible
    expect(header_texts).toContain(`Model`)
    expect(header_texts).toContain(`Prec`)
    expect(header_texts).toContain(`Acc`)
  })

  it(`filters non-compliant models`, async () => {
    let model_filter = $state(() => false) // initially show no models
    mount(MetricsTable, {
      target: document.body,
      props: { discovery_set: `unique_prototypes`, model_filter },
    })

    const initial_rows = document.body.querySelectorAll(`tbody tr`).length
    expect(initial_rows).toBe(0)

    model_filter = () => true // now show all models
    await tick()

    const rows_with_non_compliant = document.body.querySelectorAll(`tbody tr`).length
    expect(rows_with_non_compliant).toBeGreaterThan(0)
  })

  it(`opens and closes prediction files modal`, async () => {
    mount(MetricsTable, {
      target: document.body,
      props: {
        discovery_set: `unique_prototypes`,
      },
    })

    const pred_files_btn = document.body.querySelector(
      `.pred-files-btn`,
    ) as HTMLButtonElement
    const modal = document.body.querySelector(`dialog`)
    expect(modal?.open).toBe(false)

    // Open modal
    pred_files_btn?.click()
    await tick()
    expect(modal?.open).toBe(true)

    // Check modal content
    const modal_title = modal?.querySelector(`h3`)
    expect(modal_title?.textContent).toMatch(/Download prediction files for/)

    const file_links = modal?.querySelectorAll(`a`)
    expect(file_links?.length).toBeGreaterThan(1)

    // Close modal with × button
    const close_btn = modal?.querySelector(`.close-btn`) as HTMLButtonElement
    close_btn?.click()
    await tick()
    expect(modal?.open).toBe(false)

    // Check modal can be closed with Escape key
    pred_files_btn?.click()
    await tick()
    expect(modal?.open).toBe(true)

    window.dispatchEvent(new KeyboardEvent(`keydown`, { key: `Escape` }))
    await tick()
    expect(modal?.open).toBe(false)
  })

  it.each([
    {
      name: `sticky columns only`,
      col_filter: (col: HeatmapColumn) => col.sticky === true,
      expected_headers: [`Model`],
    },
    {
      name: `specific columns`,
      col_filter: (col: HeatmapColumn) => [`Model`, `F1`].includes(col.label),
      expected_headers: [`Model`, `F1`],
    },
    {
      name: `Model always first`,
      col_filter: (col: HeatmapColumn) => [`F1`, `Model`, `DAF`].includes(col.label),
      expected_headers: [`Model`, `F1`, `DAF`],
    },
  ])(`handles col_filter: $name`, async ({ col_filter, expected_headers }) => {
    mount(MetricsTable, {
      target: document.body,
      props: { discovery_set: `unique_prototypes`, col_filter },
    })

    const headers = [...document.body.querySelectorAll(`th`)]
    expect(headers.length).toBe(expected_headers.length)
    expect(headers.map((h) => h.textContent?.split(` `)[0])).toEqual(expected_headers)
  })

  it.each([
    {
      model_filter: (model: ModelData) => model.model_name.includes(`CHG`),
      col_filter: (col: HeatmapColumn) => col.label === `Model` || col.label === `F1`,
      expected_model_match: `CHG`,
      expected_headers: [`Model`, `F1`],
    },
    {
      model_filter: (model: ModelData) => model.model_name.includes(`MACE`),
      col_filter: (col: HeatmapColumn) => [`Model`, `DAF`].includes(col.label),
      expected_model_match: `MACE`,
      expected_headers: [`Model`, `DAF`],
    },
  ])(
    `combines filters: $expected_model_match models with $expected_headers`,
    async ({ model_filter, col_filter, expected_model_match, expected_headers }) => {
      mount(MetricsTable, {
        target: document.body,
        props: { discovery_set: `unique_prototypes`, model_filter, col_filter },
      })

      const headers = [...document.body.querySelectorAll(`th`)]
      expect(headers.map((h) => h.textContent?.split(` `)[0])).toEqual(expected_headers)

      const rows = document.body.querySelectorAll(`tbody tr`)
      rows.forEach((row) => {
        const model_cell = row.querySelector(`td`)
        expect(model_cell?.textContent).toContain(expected_model_match)
      })
    },
  )

  it(`updates table when col_filter changes`, async () => {
    let col_filter = $state((col: HeatmapColumn) => [`Model`, `F1`].includes(col.label))
    mount(MetricsTable, {
      target: document.body,
      props: {
        discovery_set: `unique_prototypes`,
        col_filter,
      },
    })

    // Initially should show only Model and F1
    let headers = document.body.querySelectorAll(`th`)
    expect(headers.length).toBe(2)

    // Add DAF column
    col_filter = (col) => [`Model`, `F1`, `DAF`].includes(col.label)
    await tick()

    headers = document.body.querySelectorAll(`th`)
    expect(headers.length).toBe(3)
    expect([...headers].map((h) => h.textContent?.split(` `)[0])).toEqual([
      `Model`,
      `F1`,
      `DAF`,
    ])
  })
})
