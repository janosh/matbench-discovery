import { mount, tick } from 'svelte'
import { describe, expect, it } from 'vitest'
import TableColumnToggleMenu from '../src/lib/TableColumnToggleMenu.svelte'

describe(`TableColumnToggleMenu`, () => {
  const visible_cols = {
    'Column 1': true,
    'Column 2': false,
    'Column 3': true,
  }

  it(`renders correctly with initial state`, () => {
    mount(TableColumnToggleMenu, {
      target: document.body,
      props: { visible_cols, column_panel_open: false },
    })

    // Check if the summary text is rendered
    const summary = document.body.querySelector(`summary`)
    expect(summary?.textContent).toContain(`Columns`)

    // Check if all column labels are rendered
    Object.keys(visible_cols).forEach((col) => {
      expect(document.body.textContent).toContain(col)
    })

    // Check if checkboxes reflect initial state
    const checkboxes = document.body.querySelectorAll(`input[type="checkbox"]`)
    expect(checkboxes).toHaveLength(3)
    expect((checkboxes[0] as HTMLInputElement).checked).toBe(true) // Column 1
    expect((checkboxes[1] as HTMLInputElement).checked).toBe(false) // Column 2
    expect((checkboxes[2] as HTMLInputElement).checked).toBe(true) // Column 3
  })

  it(`toggles column visibility when checkbox clicked`, async () => {
    mount(TableColumnToggleMenu, {
      target: document.body,
      props: { visible_cols, column_panel_open: false },
    })

    const checkboxes = document.body.querySelectorAll(`input[type="checkbox"]`)
    ;(checkboxes[0] as HTMLInputElement).click()
    await tick()

    // Check if the visible_cols prop was updated
    expect(visible_cols[`Column 1`]).toBe(false)
  })

  it(`opens and closes panel`, async () => {
    mount(TableColumnToggleMenu, {
      target: document.body,
      props: { visible_cols, column_panel_open: true },
    })

    const details = document.body.querySelector(`details`)
    expect(details?.open).toBe(false)

    // Click to open panel
    details?.click()
    await tick()
    expect(details?.open).toBe(true)

    // Click outside to close panel (click_outside action)
    document.body.click()
    await tick()
    expect(details?.open).toBe(false)
  })

  it(`handles HTML in column names`, () => {
    const html_cols = {
      'Column with <sub>subscript</sub>': true,
      'Column with <sup>superscript</sup>': false,
    }

    mount(TableColumnToggleMenu, {
      target: document.body,
      props: { visible_cols: html_cols, column_panel_open: false },
    })

    // Check if HTML is rendered correctly
    expect(document.body.querySelector(`sub`)).toBeTruthy()
    expect(document.body.querySelector(`sup`)).toBeTruthy()
  })
})
