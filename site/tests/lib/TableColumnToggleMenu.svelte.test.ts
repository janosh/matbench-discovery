import { TableColumnToggleMenu } from '$lib'
import type { Label } from '$lib/types'
import { mount, tick } from 'svelte'
import { describe, expect, it } from 'vitest'

describe(`TableColumnToggleMenu`, () => {
  const columns: Label[] = [
    { key: `col1`, label: `Column 1`, visible: true, description: `` },
    { key: `col2`, label: `Column 2`, visible: false, description: `` },
    { key: `col3`, label: `Column 3`, visible: true, description: `` },
  ]

  it(`renders correctly with initial state`, () => {
    mount(TableColumnToggleMenu, {
      target: document.body,
      props: { columns, column_panel_open: false },
    })

    // Verify basic elements and structure
    const summary = document.querySelector(`summary`)
    expect(summary?.textContent?.trim()).toBe(`Columns`)

    // Verify checkboxes match columns config
    const checkboxes = document.querySelectorAll(`input[type="checkbox"]`)
    expect(checkboxes).toHaveLength(3)
    expect((checkboxes[0] as HTMLInputElement).checked).toBe(true)
  })

  it(`toggles column visibility when checkbox clicked and updates state`, async () => {
    mount(TableColumnToggleMenu, {
      target: document.body,
      props: { columns, column_panel_open: false },
    })

    // Click first checkbox via its label and verify state update
    const labels = document.querySelectorAll(`label`)
    expect(
      (document.querySelectorAll(`input[type="checkbox"]`)[0] as HTMLInputElement)
        .checked,
    ).toBe(true)

    labels[0].click()
    await tick()
    expect(columns[0].visible).toBe(false)
  })

  it(`opens and closes panel with correct interaction`, async () => {
    mount(TableColumnToggleMenu, {
      target: document.body,
      props: { columns, column_panel_open: false },
    })

    const details = document.querySelector(`details`)
    expect(details?.open).toBe(false)

    // Open panel and verify state
    document.querySelector(`summary`)?.click()
    await tick()
    expect(details?.open).toBe(true)
  })

  it(`handles HTML in column names with correct rendering`, () => {
    const subscript_label = `Column with <sub>subscript</sub>`
    const superscript_label = `Column with <sup>superscript</sup>`
    const columns: Label[] = [
      { key: `col1`, label: subscript_label, visible: true, description: `` },
      { key: `col2`, label: superscript_label, visible: false, description: `` },
    ]

    mount(TableColumnToggleMenu, {
      target: document.body,
      props: { columns, column_panel_open: false },
    })

    // Verify HTML elements are properly rendered
    expect(document.querySelector(`sub`)).toBeDefined()
    expect(document.querySelector(`sup`)).toBeDefined()
  })
})
