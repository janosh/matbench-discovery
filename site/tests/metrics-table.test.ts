import MetricsTable from '$lib/MetricsTable.svelte'
import { tick } from 'svelte'
import { beforeEach, describe, expect, it } from 'vitest'

describe(`MetricsTable`, () => {
  beforeEach(() => {
    document.body.innerHTML = ``
  })

  it(`renders with default props`, () => {
    new MetricsTable({
      target: document.body,
      props: {
        discovery_set: `unique_prototypes`,
        show_non_compliant: false,
        show_energy_only: false,
        show_metadata: true,
        hide_cols: [],
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
    const component = new MetricsTable({
      target: document.body,
      props: {
        discovery_set: `unique_prototypes`,
        show_metadata: true,
      },
    })

    // Check metadata columns are visible
    let header_texts = [...document.body.querySelectorAll(`th`)].map((h) =>
      h.textContent?.trim(),
    )
    const metadata_cols = [`Training Set`, `Params`, `Targets`, `Date Added`, `Links`]
    expect(header_texts).toEqual(expect.arrayContaining(metadata_cols))

    // Hide metadata columns
    await component.$set({ show_metadata: false })
    await tick()

    // Check none of the metadata columns are hidden
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

  it(`hides specified columns`, async () => {
    new MetricsTable({
      target: document.body,
      props: {
        discovery_set: `unique_prototypes`,
        hide_cols: [`F1`, `DAF`],
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
    const component = new MetricsTable({
      target: document.body,
      props: {
        discovery_set: `unique_prototypes`,
        show_non_compliant: false,
      },
    })

    const initial_rows = document.body.querySelectorAll(`tbody tr`).length
    expect(initial_rows).toBeGreaterThan(7)

    await component.$set({ show_non_compliant: true })
    await tick()

    const rows_with_non_compliant = document.body.querySelectorAll(`tbody tr`).length
    expect(rows_with_non_compliant).toBeGreaterThan(initial_rows)
  })

  it(`opens and closes prediction files modal`, async () => {
    new MetricsTable({
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
})
