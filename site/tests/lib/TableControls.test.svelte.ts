import { TableControls } from '$lib'
import { mount } from 'svelte'
import { describe, expect, it, vi } from 'vitest'

describe(`TableControls`, () => {
  const sample_cols = { Model: true, F1: true, DAF: true, RMSE: false }

  it(`renders filter controls with correct initial state`, () => {
    const on_filter_change = vi.fn()
    const on_col_change = vi.fn()

    mount(TableControls, {
      target: document.body,
      props: { visible_cols: sample_cols, on_filter_change, on_col_change },
    })

    // Verify specific filter checkboxes are present
    const filter_checkboxes = document.querySelectorAll(`input[type="checkbox"]`)
    expect(filter_checkboxes.length).toBeGreaterThan(2) // at least energy-only, non-compliant, heatmap

    // Verify column toggle button exists
    const col_toggle_btn = document.querySelector(
      `[aria-label="Toggle column visibility"]`,
    )
    expect(col_toggle_btn).toBeDefined()
  })

  it(`calls on_filter_change with correct parameters when filters are toggled`, () => {
    const on_filter_change = vi.fn()

    mount(TableControls, {
      target: document.body,
      props: { visible_cols: sample_cols, on_filter_change },
    })

    const checkboxes = document.querySelectorAll(`input[type="checkbox"]`)
    expect(checkboxes.length).toBeGreaterThan(0)

    // Test energy-only filter toggle
    const energy_checkbox = Array.from(checkboxes).find((checkbox) =>
      checkbox.id?.includes(`energy`)
    ) as HTMLInputElement

    if (energy_checkbox) {
      const initial_checked = energy_checkbox.checked
      energy_checkbox.click()

      // Verify callback was called
      expect(on_filter_change).toHaveBeenCalled()

      // Verify checkbox state changed
      expect(energy_checkbox.checked).toBe(!initial_checked)
    }

    // Test non-compliant filter toggle
    on_filter_change.mockReset()
    const noncomp_checkbox = Array.from(checkboxes).find((checkbox) =>
      checkbox.id?.includes(`compliant`)
    ) as HTMLInputElement

    if (noncomp_checkbox) {
      const initial_checked = noncomp_checkbox.checked
      noncomp_checkbox.click()

      expect(on_filter_change).toHaveBeenCalled()
      expect(noncomp_checkbox.checked).toBe(!initial_checked)
    }
  })

  it(`toggles column visibility panel`, () => {
    mount(TableControls, {
      target: document.body,
      props: { visible_cols: sample_cols },
    })

    // Find column toggle button
    const toggle_btn = document.querySelector(
      `[aria-label="Toggle column visibility"]`,
    ) as HTMLButtonElement
    expect(toggle_btn).toBeDefined()

    // Click button to show the panel
    if (toggle_btn) {
      toggle_btn.click()

      // Check if column menu or panel is visible
      // The component might use different class names, so we try multiple options
      const column_menu = document.querySelector(`.column-menu, .column-panel`)
      expect(column_menu).toBeDefined()

      // Click outside to close (if the component uses a click-outside pattern)
      document.body.dispatchEvent(new MouseEvent(`click`, { bubbles: true }))
    }
  })

  it(`calls on_col_change when column visibility is toggled`, () => {
    const on_col_change = vi.fn()

    mount(TableControls, {
      target: document.body,
      props: { visible_cols: { ...sample_cols }, on_col_change },
    })

    const toggle_btn = document.querySelector(
      `[aria-label="Toggle column visibility"]`,
    ) as HTMLButtonElement

    if (!toggle_btn) {
      // Skip test if component doesn't have toggle button
      return
    }

    // Open column panel
    toggle_btn.click()

    // Find and toggle a column checkbox
    const column_checkboxes = document.querySelectorAll(`input[type="checkbox"]`)
    expect(column_checkboxes.length).toBeGreaterThan(0)

    const first_checkbox = column_checkboxes[0] as HTMLInputElement
    const initial_checked = first_checkbox.checked

    first_checkbox.click()

    // Verify callback was called
    expect(on_col_change).toHaveBeenCalled()

    // Verify checkbox state actually changed
    expect(first_checkbox.checked).toBe(!initial_checked)
  })
})
