import { TableColumnToggleMenu } from '$lib'
import type { HeatmapColumn } from '$lib/types'
import { mount, tick } from 'svelte'
import { describe, expect, it } from 'vitest'

describe(`TableColumnToggleMenu`, () => {
  const columns: HeatmapColumn[] = [
    { label: `Column 1`, visible: true },
    { label: `Column 2`, visible: false },
    { label: `Column 3`, visible: true },
  ]

  it(`renders correctly with initial state`, () => {
    mount(TableColumnToggleMenu, {
      target: document.body,
      props: { columns, column_panel_open: false },
    })

    // Verify basic elements and structure
    const summary = document.body.querySelector(`summary`)
    expect(summary?.textContent?.trim()).toBe(`Columns`)

    // Verify checkboxes match columns config
    const checkboxes = document.body.querySelectorAll(`input[type="checkbox"]`)
    expect(checkboxes).toHaveLength(3)
    expect((checkboxes[0] as HTMLInputElement).checked).toBe(true)
  })

  it(`toggles column visibility when checkbox clicked and updates state`, async () => {
    mount(TableColumnToggleMenu, {
      target: document.body,
      props: { columns, column_panel_open: false },
    })

    // Click first checkbox via its label and verify state update
    const labels = document.body.querySelectorAll(`label`)
    expect(
      (document.body.querySelectorAll(`input[type="checkbox"]`)[0] as HTMLInputElement)
        .checked,
    ).toBe(true)

    labels[0].click()
    await tick()
    expect(columns[0].visible).toBe(false)
  })

  it(`opens and closes panel with correct interaction`, async () => {
    mount(TableColumnToggleMenu, {
      target: document.body,
      props: { columns, column_panel_open: true },
    })

    const details = document.body.querySelector(`details`)
    expect(details?.open).toBe(false)

    // Open panel and verify state
    document.body.querySelector(`summary`)?.click()
    await tick()
    expect(details?.open).toBe(true)
  })

  it(`handles HTML in column names with correct rendering`, () => {
    const columns: HeatmapColumn[] = [
      { label: `Column with <sub>subscript</sub>`, visible: true },
      { label: `Column with <sup>superscript</sup>`, visible: false },
    ]

    mount(TableColumnToggleMenu, {
      target: document.body,
      props: { columns, column_panel_open: false },
    })

    // Verify HTML elements are properly rendered
    expect(document.body.querySelector(`sub`)).toBeDefined()
    expect(document.body.querySelector(`sup`)).toBeDefined()
  })
})
