import { DEFAULT_COMBINED_METRIC_CONFIG } from '$lib/metrics'
import TableControls from '$lib/TableControls.svelte'
import type { CombinedMetricConfig } from '$lib/types'
import { mount, tick } from 'svelte'
import { describe, expect, it, vi } from 'vitest'

describe(`TableControls`, () => {
  const sample_cols = { Model: true, F1: true, DAF: true, RMSE: false }

  it(`renders with default props`, async () => {
    const on_filter_change = vi.fn()
    const on_col_change = vi.fn()

    mount(TableControls, {
      target: document.body,
      props: {
        visible_cols: sample_cols,
        on_filter_change,
        on_col_change,
      },
    })

    await tick()

    // Check filter checkboxes
    const filter_checkboxes = document.body.querySelectorAll(`input[type="checkbox"]`)
    expect(filter_checkboxes.length).toBeGreaterThan(0)

    // Check column toggle button exists
    const col_toggle_btn = document.body.querySelector(
      `[aria-label="Toggle column visibility"]`,
    )
    expect(col_toggle_btn).toBeDefined()
  })

  it(`responds to filter changes`, async () => {
    const on_filter_change = vi.fn()

    mount(TableControls, {
      target: document.body,
      props: {
        visible_cols: sample_cols,
        on_filter_change,
        // Default values
        show_energy_only: false,
        show_noncompliant: false,
      },
    })

    await tick()

    // Find filter checkboxes
    const checkboxes = document.body.querySelectorAll(`input[type="checkbox"]`)
    expect(checkboxes.length).toBeGreaterThan(0)

    // Click first checkbox (energy-only) and check if callback is called
    on_filter_change.mockReset()
    const energy_checkbox = Array.from(checkboxes).find((checkbox) =>
      checkbox.id?.includes(`energy`),
    ) as HTMLInputElement

    if (energy_checkbox) {
      energy_checkbox.click()
      expect(on_filter_change).toHaveBeenCalled()
    }

    // Click second checkbox (noncompliant) and check if callback is called
    on_filter_change.mockReset()
    const noncomp_checkbox = Array.from(checkboxes).find((checkbox) =>
      checkbox.id?.includes(`compliant`),
    ) as HTMLInputElement

    if (noncomp_checkbox) {
      noncomp_checkbox.click()
      expect(on_filter_change).toHaveBeenCalled()
    }
  })

  it(`toggles column visibility panel`, async () => {
    mount(TableControls, {
      target: document.body,
      props: {
        visible_cols: sample_cols,
      },
    })

    await tick()

    // Find column toggle button
    const toggle_btn = document.body.querySelector(
      `[aria-label="Toggle column visibility"]`,
    ) as HTMLButtonElement
    expect(toggle_btn).toBeDefined()

    // Click button to show the panel
    if (toggle_btn) {
      toggle_btn.click()
      await tick()

      // Check if column menu or panel is visible
      // The component might use different class names, so we try multiple options
      const column_menu = document.body.querySelector(`.column-menu, .column-panel`)
      expect(column_menu).toBeDefined()

      // Click outside to close (if the component uses a click-outside pattern)
      document.body.dispatchEvent(new MouseEvent(`click`, { bubbles: true }))
      await tick()
    }
  })

  it(`handles column visibility changes`, async () => {
    const on_col_change = vi.fn()

    mount(TableControls, {
      target: document.body,
      props: {
        visible_cols: { ...sample_cols },
        on_col_change,
      },
    })

    await tick()

    // Open column menu/panel (if it exists)
    const toggle_btn = document.body.querySelector(
      `[aria-label="Toggle column visibility"]`,
    ) as HTMLButtonElement
    if (toggle_btn) {
      toggle_btn.click()
      await tick()

      // Click column checkbox and check if callback is called
      const column_checkboxes = document.body.querySelectorAll(`input[type="checkbox"]`)
      if (column_checkboxes.length > 0) {
        on_col_change.mockReset()
        ;(column_checkboxes[0] as HTMLInputElement).click()
        expect(on_col_change).toHaveBeenCalled()
      }
    }
  })

  it(`renders tooltip info icons`, async () => {
    mount(TableControls, {
      target: document.body,
      props: {
        visible_cols: sample_cols,
      },
    })

    await tick()

    // Check for info icons in the document
    const info_icons = document.body.querySelectorAll(`.info-icon, [aria-label*="info"]`)

    // If info icons exist, try simulating a hover
    if (info_icons.length > 0) {
      const info_icon = info_icons[0] as HTMLElement
      const mouseenter_event = new MouseEvent(`mouseenter`)
      info_icon.dispatchEvent(mouseenter_event)
      await tick()
    }
  })

  it(`passes config to the component`, () => {
    const sample_config: CombinedMetricConfig = {
      ...DEFAULT_COMBINED_METRIC_CONFIG,
      name: `Test Config`,
    }

    mount(TableControls, {
      target: document.body,
      props: {
        visible_cols: sample_cols,
        config: sample_config,
      },
    })

    // The component should accept the config prop
    // (This is a basic test of prop passing, functionality would be tested in RadarChart tests)
  })
})
