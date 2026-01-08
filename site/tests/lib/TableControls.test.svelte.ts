import { type Label, TableControls } from '$lib'
import { mount } from 'svelte'
import { describe, expect, it, vi } from 'vitest'
import { doc_query } from '../index.ts'

describe(`TableControls`, () => {
  const sample_columns: Label[] = [
    { key: `model`, label: `Model`, description: `Model name`, visible: true },
    { key: `f1`, label: `F1`, description: `F1 Score`, visible: true },
    { key: `daf`, label: `DAF`, description: `DAF Score`, visible: true },
    { key: `rmse`, label: `RMSE`, description: `RMSE`, visible: false },
  ]

  // Helper to find checkbox by parent label text
  const find_checkbox_by_label = (text: string): HTMLInputElement | null => {
    const labels = document.querySelectorAll(`label`)
    const label = Array.from(labels).find((lbl) => lbl.textContent?.includes(text))
    return doc_query<HTMLInputElement>(`input[type="checkbox"]`, label)
  }

  it(`renders filter controls with correct initial state`, () => {
    mount(TableControls, {
      target: document.body,
      props: { on_filter_change: vi.fn() },
    })

    // Verify filter checkboxes are present
    expect(document.querySelectorAll(`input[type="checkbox"]`).length).toBeGreaterThan(2)

    // Verify column toggle exists
    expect(document.querySelector(`.column-toggles summary`)).toBeTruthy()
  })

  it(`calls on_filter_change when energy-only filter is toggled`, () => {
    const on_filter_change = vi.fn()

    mount(TableControls, {
      target: document.body,
      props: { on_filter_change, show_energy_only_toggle: true },
    })

    const energy_checkbox = find_checkbox_by_label(`Energy-only`)
    expect(energy_checkbox).toBeTruthy()
    if (!energy_checkbox) return

    const initial_checked = energy_checkbox.checked
    energy_checkbox.click()

    expect(on_filter_change).toHaveBeenCalledWith(!initial_checked, false)
    expect(energy_checkbox.checked).toBe(!initial_checked)
  })

  it(`toggles compliance filter checkboxes`, () => {
    mount(TableControls, { target: document.body })

    const compliant_checkbox = find_checkbox_by_label(`Compliant`)
    const noncompliant_checkbox = find_checkbox_by_label(`Non-compliant`)

    expect(compliant_checkbox).toBeTruthy()
    expect(noncompliant_checkbox).toBeTruthy()
    if (!compliant_checkbox || !noncompliant_checkbox) return

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

  it(`opens and closes column visibility panel`, () => {
    mount(TableControls, {
      target: document.body,
      props: { columns: sample_columns },
    })

    const toggle_btn = doc_query(`.column-toggles summary`)
    const details = doc_query<HTMLDetailsElement>(`.column-toggles`)
    expect(toggle_btn).toBeTruthy()
    expect(details).toBeTruthy()

    // Should start closed
    expect(details.open).toBe(false)

    // Open panel
    toggle_btn.click()
    expect(details.open).toBe(true)

    // Verify column menu is visible
    expect(document.querySelector(`.column-menu`)).toBeTruthy()

    // Close panel
    toggle_btn.click()
    expect(details.open).toBe(false)
  })

  it(`toggles column visibility checkboxes`, () => {
    mount(TableControls, {
      target: document.body,
      props: { columns: [...sample_columns] },
    })

    const toggle_btn = doc_query(`.column-toggles summary`)
    toggle_btn.click()

    const column_menu = document.querySelector(`.column-menu`)
    expect(column_menu).toBeTruthy()

    const column_checkboxes = column_menu?.querySelectorAll(`input[type="checkbox"]`)
    expect(column_checkboxes?.length).toBe(sample_columns.length)

    // Test toggling first checkbox (Model column, initially visible)
    const first_checkbox = column_checkboxes?.[0] as HTMLInputElement
    expect(first_checkbox?.checked).toBe(true)

    first_checkbox.click()
    expect(first_checkbox.checked).toBe(false)

    first_checkbox.click()
    // Verify checkbox state changed
    expect(first_checkbox.checked).toBe(true)
  })
})
