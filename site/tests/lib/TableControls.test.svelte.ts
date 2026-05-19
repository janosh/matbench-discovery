import { type Label, TableControls } from '$lib'
import { mount, tick } from 'svelte'
import { describe, expect, it, vi } from 'vitest'
import { doc_query } from '../index'

describe(`TableControls`, () => {
  const sample_columns: Label[] = [
    { key: `model`, label: `Model`, description: `Model name`, visible: true },
    { key: `f1`, label: `F1`, description: `F1 Score`, visible: true },
    { key: `daf`, label: `DAF`, description: `DAF Score`, visible: true },
    { key: `rmse`, label: `RMSE`, description: `RMSE`, visible: false },
  ]

  const find_checkbox_by_label = (text: string): HTMLInputElement => {
    const labels = document.querySelectorAll(`label`)
    const label = [...labels].find((lbl) => lbl.textContent?.includes(text))
    if (!label) throw new Error(`No checkbox label found containing: ${text}`)
    return doc_query<HTMLInputElement>(`input[type="checkbox"]`, label)
  }

  it(`renders filter controls with correct initial state`, () => {
    mount(TableControls, {
      target: document.body,
      props: { on_filter_change: vi.fn() },
    })

    // Verify filter checkboxes are present
    const labels = [...document.querySelectorAll(`label`)].map((label) =>
      label.textContent?.replace(/\s+/g, ` `).trim(),
    )
    expect(labels).toContain(`Compliant models`)
    expect(labels).toContain(`Non-compliant models`)
    expect(labels).toContain(`Heatmap`)
  })

  it(`calls on_filter_change when energy-only filter is toggled`, () => {
    const on_filter_change = vi.fn()

    mount(TableControls, {
      target: document.body,
      props: { on_filter_change, show_energy_only_toggle: true },
    })

    const energy_checkbox = find_checkbox_by_label(`Energy-only`)

    const initial_checked = energy_checkbox.checked
    energy_checkbox.click()

    expect(on_filter_change).toHaveBeenCalledExactlyOnceWith(!initial_checked, false)
    expect(energy_checkbox.checked).toBe(!initial_checked)
  })

  it(`toggles compliance filter checkboxes`, () => {
    mount(TableControls, { target: document.body })

    const compliant_checkbox = find_checkbox_by_label(`Compliant`)
    const noncompliant_checkbox = find_checkbox_by_label(`Non-compliant`)

    // Both should start checked
    expect(compliant_checkbox.checked).toBe(true)
    expect(noncompliant_checkbox.checked).toBe(true)

    // Toggle non-compliant off
    noncompliant_checkbox.click()
    expect(noncompliant_checkbox.checked).toBe(false)

    // Toggle back on
    noncompliant_checkbox.click()
    expect(noncompliant_checkbox.checked).toBe(true)
  })

  it(`opens and closes column visibility panel`, async () => {
    mount(TableControls, {
      target: document.body,
      props: { columns: sample_columns },
    })
    await tick()

    const toggle_btn = doc_query(`.column-toggles summary`)
    const details = doc_query<HTMLDetailsElement>(`.column-toggles`)
    expect(toggle_btn.textContent?.trim()).toBe(`Columns`)

    // Should start closed
    expect(details.open).toBe(false)

    // Open panel
    toggle_btn.click()
    expect(details.open).toBe(true)

    // Verify column menu is visible
    expect(doc_query(`.column-menu`).querySelectorAll(`input[type="checkbox"]`)).toHaveLength(
      sample_columns.length,
    )

    // Close panel
    toggle_btn.click()
    expect(details.open).toBe(false)
  })

  it(`toggles column visibility checkboxes`, async () => {
    mount(TableControls, {
      target: document.body,
      props: { columns: [...sample_columns] },
    })
    await tick()

    const toggle_btn = doc_query(`.column-toggles summary`)
    toggle_btn.click()

    const column_menu = doc_query(`.column-menu`)
    expect(column_menu.getAttribute(`role`)).toBe(`group`)

    const column_checkboxes = column_menu.querySelectorAll<HTMLInputElement>(`input[type="checkbox"]`)
    expect(column_checkboxes).toHaveLength(sample_columns.length)
    const checkbox_labels = [...column_menu.querySelectorAll(`label`)].map((label) =>
      label.textContent?.trim(),
    )
    expect(checkbox_labels).toStrictEqual(sample_columns.map((column) => column.label))

    // Test toggling first checkbox (Model column, initially visible)
    const [first_checkbox] = column_checkboxes
    expect(first_checkbox.checked).toBe(true)

    first_checkbox.click()
    expect(first_checkbox.checked).toBe(false)

    first_checkbox.click()
    // Verify checkbox state changed
    expect(first_checkbox.checked).toBe(true)
  })
})
