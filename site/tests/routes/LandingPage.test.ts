import Page from '$routes/+page.svelte'
import { mount, tick } from 'svelte'
import { beforeEach, describe, expect, it } from 'vitest'
import { doc_query } from '../index'

describe(`Landing Page`, () => {
  beforeEach(() => {
    mount(Page, { target: document.body })
  })

  it(`renders discovery set toggle buttons`, () => {
    const buttons = document.querySelectorAll(`.selection-toggle button`)
    expect(buttons).toHaveLength(3) // 3 from test set select

    const button_texts = [...buttons].map((btn) => btn.textContent?.trim())
    expect(button_texts).toStrictEqual([
      `Full Test Set`,
      `Unique Prototypes`,
      `10k Most Stable`,
    ])
  })

  it(`toggles discovery set when clicking buttons`, async () => {
    const buttons = [...document.querySelectorAll<HTMLButtonElement>(`.selection-toggle button`)]
    const [full_test_btn, unique_protos_btn] = buttons

    // Initially Unique Prototypes should be active
    expect(unique_protos_btn.classList.contains(`active`)).toBe(true)
    expect(full_test_btn.classList.contains(`active`)).toBe(false)

    // Click Full Test Set button
    full_test_btn.click()
    await tick()

    expect(unique_protos_btn.classList.contains(`active`)).toBe(false)
    expect(full_test_btn.classList.contains(`active`)).toBe(true)
  })

  it(`toggles column visibility panel`, () => {
    const columns_btn = doc_query<HTMLButtonElement>(`details.column-toggles summary`)
    const column_menu = doc_query(`.column-menu`)
    const details = doc_query<HTMLDetailsElement>(`details.column-toggles`)

    // Column menu should be hidden initially
    expect(column_menu.closest(`details`)).toBe(details)
    expect(details.open).toBe(false)

    // Click columns button to open menu
    columns_btn.click()
    expect(details.open).toBe(true)

    // Click outside to close menu
    columns_btn.click()
    expect(details.open).toBe(false)
  })

  it(`toggles non-compliant models`, async () => {
    const toggle = doc_query<HTMLInputElement>(
      `.table-controls label input[type="checkbox"]`,
    )
    expect(toggle.parentElement?.textContent).toContain(`Compliant models`)

    // Should be unchecked by default
    expect(toggle.checked).toBe(true)
    // get number of table rows
    const n_models_on_load = document.querySelectorAll(`tbody tr`).length

    // Click to show non-compliant models
    toggle.click()
    expect(toggle.checked).toBe(false)
    await tick()
    const n_all_models = document.querySelectorAll(`tbody tr`).length
    expect(n_all_models).toBeLessThan(n_models_on_load)
  })

  it(`updates column visibility when toggling checkboxes`, async () => {
    const columns_btn = doc_query<HTMLButtonElement>(`details.column-toggles summary`)
    columns_btn.click()
    // Table should reflect column visibility changes
    let f1_cells = document.querySelectorAll(`th, td`)
    let has_f1_column = [...f1_cells].some((cell) =>
      cell.textContent?.includes(`F1`),
    )
    expect(has_f1_column).toBe(true)

    const checkboxes = document.querySelectorAll<HTMLInputElement>(
      `.column-menu input[type="checkbox"]`,
    )
    const f1_checkbox = [...checkboxes].find((cb) =>
      cb.parentElement?.textContent?.includes(`F1`),
    )
    if (!f1_checkbox) throw new Error(`F1 checkbox not found`)
    expect(f1_checkbox.checked).toBe(true)

    // Uncheck F1
    f1_checkbox.click()
    await tick()
    expect(f1_checkbox.checked).toBe(false)

    // Table should reflect column visibility changes
    f1_cells = document.querySelectorAll(`th, td`)
    has_f1_column = [...f1_cells].some((cell) => cell.textContent?.includes(`F1`))
    expect(has_f1_column).toBe(false)
  })

  it(`displays best model information`, () => {
    const best_model_info = doc_query(`#best-report`)
    expect(best_model_info.textContent).toMatch(/highest F1 score/)
    expect(best_model_info.textContent).toMatch(/discovery acceleration factor/)
  })

  it(`renders table downloads section`, () => {
    // Check that the download buttons section exists
    const download_section = doc_query(`.downloads`)

    // Check that it contains SVG, PNG, CSV, Excel buttons plus RSS link
    const download_buttons = download_section.querySelectorAll(`.download-btn`)
    expect(download_buttons).toHaveLength(5)

    const buttons = [...download_buttons].map((btn) =>
      btn.textContent?.trim(),
    )
    expect(buttons).toStrictEqual([`SVG`, `PNG`, `CSV`, `Excel`, `RSS`])
  })

  it(`displays valid metric values`, () => {
    const text = doc_query(`#best-report`).textContent ?? ``

    // Extract F1 and DAF values using regex
    const f1_match = /F1 score of ([\d.]+)/.exec(text)
    const daf_match = /DAF\) of ([\d.]+)/.exec(text)

    if (f1_match && daf_match) {
      const f1_value = parseFloat(f1_match[1])
      const daf_value = parseFloat(daf_match[1])

      expect(f1_value > 0 && f1_value < 1, `F1=${f1_value} is out of range`).toBe(true)
      expect(daf_value > 0 && daf_value < 10, `DAF=${daf_value} is out of range`).toBe(
        true,
      )
    } else {
      throw new Error(`Could not find F1 or DAF values in text: ${text}`)
    }
  })
})
