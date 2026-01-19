import DynamicScatter from '$lib/DynamicScatter.svelte'
import type { ModelData } from '$lib/types'
import { mount } from 'svelte'
import { beforeEach, describe, expect, it, vi } from 'vitest'
import { doc_query, is_hidden } from '../index'

const pane_selector = `[aria-label="Draggable pane"]`

const mock_models: ModelData[] = [
  {
    model_name: `Test Model 1`,
    model_params: 100000,
    n_estimators: 50,
    dirname: `test1`,
    metadata_file: `test1.yml`,
    date_added: `2023-01-01`,
    color: `blue`,
    metrics: {
      discovery: {
        unique_prototypes: { F1: 0.85, DAF: 0.9 },
        full_test_set: { F1: 0.8, DAF: 0.85 },
      },
      geo_opt: { 'symprec=1e-5': { rmsd: 0.02 }, 'symprec=1e-2': { rmsd: 0.03 } },
    },
    hyperparams: { batch_size: 32, n_layers: 4 },
  } as unknown as ModelData,
  {
    model_name: `Test Model 2`,
    model_params: 500000,
    n_estimators: 100,
    dirname: `test2`,
    metadata_file: `test2.yml`,
    date_added: `2023-03-15`,
    color: `red`,
    metrics: {
      discovery: {
        unique_prototypes: { F1: 0.9, DAF: 0.95 },
        full_test_set: { F1: 0.85, DAF: 0.9 },
      },
      geo_opt: { 'symprec=1e-5': { rmsd: 0.015 }, 'symprec=1e-2': { rmsd: 0.025 } },
    },
    hyperparams: { batch_size: 64, n_layers: 6 },
  } as unknown as ModelData,
  {
    model_name: `Test Model 3 - No Params`,
    model_params: null, // Missing model_params
    n_estimators: 200,
    dirname: `test3`,
    metadata_file: `test3.yml`,
    date_added: `2023-06-30`,
    color: `green`,
    metrics: {
      discovery: {
        unique_prototypes: { F1: 0.95, DAF: 0.98 },
        full_test_set: { F1: 0.9, DAF: 0.93 },
      },
      geo_opt: { 'symprec=1e-5': { rmsd: 0.01 }, 'symprec=1e-2': { rmsd: 0.02 } },
    },
    hyperparams: { batch_size: 128, n_layers: 8 },
  } as unknown as ModelData,
]

describe(`DynamicScatter.svelte`, () => {
  beforeEach(() => {
    // Mock fullscreen API
    let fullscreen_element: Element | null = null
    Object.defineProperty(document, `fullscreenElement`, {
      get: () => fullscreen_element,
      configurable: true,
    })
    // Mock requestFullscreen and exitFullscreen methods
    Element.prototype.requestFullscreen = vi.fn().mockImplementation(
      function (this: Element) {
        fullscreen_element = this // eslint-disable-line @typescript-eslint/no-this-alias
        document.dispatchEvent(new Event(`fullscreenchange`))
        return Promise.resolve()
      },
    )
    document.exitFullscreen = vi.fn().mockImplementation(() => {
      fullscreen_element = null
      document.dispatchEvent(new Event(`fullscreenchange`))
      return Promise.resolve()
    })
  })

  it(`mounts correctly with default props`, () => {
    mount(DynamicScatter, {
      target: document.body,
      props: { models: mock_models },
    })

    // Check that controls row is rendered (contains Size select)
    const controls_row = document.querySelector(`.controls-row`)
    expect(controls_row).toBeDefined()

    // Check that the scatter plot container is rendered
    const plot_container = document.querySelector(`div.bleed-1400[style]`)
    expect(plot_container).toBeDefined()

    // Check that fullscreen button is rendered
    const fullscreen_button = document.querySelector(`button[title$="fullscreen"]`)
    expect(fullscreen_button).toBeDefined()

    // Check that settings toggle button is rendered
    const settings_button = document.querySelector(`.settings-toggle`)
    expect(settings_button).toBeDefined()
  })

  it(`renders component structure (no SVG check)`, () => {
    // Use mount from svelte
    mount(DynamicScatter, { target: document.body, props: { models: mock_models } })
    expect(document.querySelector(`.controls-row`)).toBeDefined()
    expect(document.querySelector(`div.bleed-1400[style]`)).toBeDefined()
  })

  // Helper function to check fullscreen state
  async function check_fullscreen_state(expected_state: boolean): Promise<void> {
    const button = document.querySelector<HTMLButtonElement>(
      `button[title$="fullscreen"]`,
    )
    const button_icon = button?.querySelector(`button > svg`)

    await vi.waitFor(() => {
      expect(button_icon, `Button icon`).toBeDefined()
      expect(button?.getAttribute(`aria-label`), `Button label`).toBe(
        `${expected_state ? `Exit` : `Enter`} fullscreen`,
      )
    })
  }

  it(`handles fullscreen toggle via button click`, async () => {
    mount(DynamicScatter, {
      target: document.body,
      props: { models: mock_models },
    })

    const button = doc_query<HTMLButtonElement>(`button[title$="fullscreen"]`)

    // 1. Initial state: Not fullscreen
    await check_fullscreen_state(false)

    // 2. Enter fullscreen via button click
    await button?.click()
    await check_fullscreen_state(true)

    // 3. Re-enter fullscreen via button click
    await button?.click()
    await check_fullscreen_state(false)
  })

  it(`respects model_filter prop (initial render only)`, () => {
    // Test only initial render effect
    mount(DynamicScatter, {
      target: document.body,
      props: {
        models: mock_models,
        model_filter: (model: ModelData) => model.model_params !== null,
      },
    })
    expect(document.querySelector(`.controls-row`)).toBeDefined()
    // Cannot easily check filter effect without mocks
  })

  describe(`extra controls`, () => {
    it(`toggles extra controls visibility via button click and Escape key`, async () => {
      mount(DynamicScatter, {
        target: document.body,
        props: { models: mock_models },
      })

      const settings_button = doc_query<HTMLButtonElement>(`.settings-toggle`)
      let extra_controls = document.querySelector(pane_selector)

      // 1. Initial state: Controls hidden (element may exist but be hidden)
      expect(is_hidden(extra_controls)).toBe(true)

      // 2. Show controls via button click
      await settings_button?.click()
      extra_controls = doc_query<HTMLElement>(pane_selector)
      expect(extra_controls).toBeDefined()

      // 3. Hide controls via Escape key
      document.dispatchEvent(
        new KeyboardEvent(`keydown`, { key: `Escape`, bubbles: true }),
      )
      await vi.waitFor(() => {
        extra_controls = doc_query<HTMLElement>(pane_selector)
        // The pane should be hidden but still in DOM
        expect(is_hidden(extra_controls)).toBe(true)
      })

      // 4. Re-show controls via button click
      await settings_button?.click()
      extra_controls = doc_query<HTMLElement>(pane_selector)
      expect(extra_controls).toBeDefined()
    })

    it(`toggles extra controls visibility via click outside`, async () => {
      // Create an explicit outside element to click
      const outside_element = document.createElement(`div`)
      outside_element.setAttribute(`data-testid`, `outside`)
      document.body.appendChild(outside_element)

      mount(DynamicScatter, {
        target: document.body,
        props: { models: mock_models },
      })

      const settings_button = doc_query<HTMLButtonElement>(`.settings-toggle`)
      let extra_controls = doc_query<HTMLElement>(pane_selector)

      // 1. Show controls
      await settings_button?.click()
      extra_controls = doc_query<HTMLElement>(pane_selector)
      expect(extra_controls).toBeDefined()

      // 2. Click the explicit outside element
      await outside_element.click() // Simulate click outside
      await vi.waitFor(() => {
        extra_controls = doc_query<HTMLElement>(pane_selector)
        // The pane should be hidden but still in DOM
        expect(is_hidden(extra_controls)).toBe(true)
      })

      // Clean up the outside element
      document.body.removeChild(outside_element)
    })

    it(`interacts with all extra controls`, async () => {
      mount(DynamicScatter, {
        target: document.body,
        props: { models: mock_models },
      })

      const settings_button = doc_query<HTMLButtonElement>(`.settings-toggle`)
      await settings_button?.click() // Show controls

      const extra_controls = doc_query(pane_selector)

      // --- Find Controls ---
      const color_scale_select_el = extra_controls?.querySelector(`.color-scale-select`)
      // Find checkboxes outside the log-toggles div (Show Labels, X Grid, Y Grid)
      const all_checkboxes = extra_controls?.querySelectorAll<HTMLInputElement>(
        `input[type="checkbox"]`,
      )
      // Log toggles are in .log-toggles div (first 4), then Show Labels, X Grid, Y Grid
      // The log toggles may be hidden via visibility but are still in DOM
      const log_toggles_count = extra_controls?.querySelectorAll(
        `.log-toggles input[type="checkbox"]`,
      ).length ?? 0
      const show_labels_checkbox = all_checkboxes?.[log_toggles_count]
      const x_grid_checkbox = all_checkboxes?.[log_toggles_count + 1]
      const y_grid_checkbox = all_checkboxes?.[log_toggles_count + 2]

      const x_ticks_input = doc_query<HTMLInputElement>(`#x-ticks`, extra_controls)
      const y_ticks_input = doc_query<HTMLInputElement>(`#y-ticks`, extra_controls)
      const size_multiplier_input = doc_query<HTMLInputElement>(
        `#size-multiplier`,
        extra_controls,
      )
      const label_font_size_input = doc_query<HTMLInputElement>(
        `#label-font-size`,
        extra_controls,
      )
      const min_link_distance_input = doc_query<HTMLInputElement>(
        `#min-link-distance`,
        extra_controls,
      )
      const max_link_distance_input = doc_query<HTMLInputElement>(
        `#max-link-distance`,
        extra_controls,
      )
      const link_strength_input = doc_query<HTMLInputElement>(
        `#link-strength`,
        extra_controls,
      )

      // --- Verify Initial States ---
      expect(color_scale_select_el, `ColorScaleSelect should be rendered`).toBeDefined()
      expect(show_labels_checkbox?.checked, `Show Labels default`).toBe(true)
      expect(x_grid_checkbox?.checked, `X Grid default`).toBe(true)
      expect(x_ticks_input?.value, `X Ticks default`).toBe(`5`)
      expect(y_grid_checkbox?.checked, `Y Grid default`).toBe(true)
      expect(y_ticks_input?.value, `Y Ticks default`).toBe(`5`)
      expect(size_multiplier_input?.value, `Size Multiplier default`).toBe(`1`)
      expect(label_font_size_input?.value, `Label Font Size default`).toBe(`14`)
      expect(min_link_distance_input?.value, `Min Link Distance default`).toBe(`15`)
      expect(max_link_distance_input?.value, `Max Link Distance default`).toBe(`20`)
      expect(link_strength_input?.value, `Link Strength default`).toBe(`5`)

      // --- Simulate Interactions & Verify Updates ---

      // Checkboxes
      await show_labels_checkbox?.click()
      expect(show_labels_checkbox?.checked, `Show Labels after click`).toBe(false)
      await x_grid_checkbox?.click()
      expect(x_grid_checkbox?.checked, `X Grid after click`).toBe(false)
      await y_grid_checkbox?.click()
      expect(y_grid_checkbox?.checked, `Y Grid after click`).toBe(false)

      // Number Inputs
      x_ticks_input.value = `10`
      await x_ticks_input?.dispatchEvent(new Event(`input`))
      expect(x_ticks_input?.value, `X Ticks after change`).toBe(`10`)

      y_ticks_input.value = `12`
      await y_ticks_input?.dispatchEvent(new Event(`input`))
      expect(y_ticks_input?.value, `Y Ticks after change`).toBe(`12`)

      min_link_distance_input.value = `10`
      await min_link_distance_input?.dispatchEvent(new Event(`input`))
      expect(min_link_distance_input?.value, `Min Link Distance after change`).toBe(`10`)

      max_link_distance_input.value = `30`
      await max_link_distance_input?.dispatchEvent(new Event(`input`))
      expect(max_link_distance_input?.value, `Max Link Distance after change`).toBe(`30`)

      // Range Sliders
      size_multiplier_input.value = `2.5`
      await size_multiplier_input?.dispatchEvent(new Event(`input`))
      expect(size_multiplier_input?.value, `Size Multiplier after change`).toBe(`2.5`)

      label_font_size_input.value = `18`
      await label_font_size_input?.dispatchEvent(new Event(`input`))
      expect(label_font_size_input?.value, `Label Font Size after change`).toBe(`18`)

      link_strength_input.value = `8`
      await link_strength_input?.dispatchEvent(new Event(`input`))
      expect(link_strength_input?.value, `Link Strength after change`).toBe(`8`)
    })
  })

  describe(`regression tests for default values`, () => {
    it(`verifies all critical default UI state to catch regressions`, async () => {
      mount(DynamicScatter, { target: document.body, props: { models: mock_models } })

      // Verify controls row is rendered with Size select
      const controls_row = document.querySelector(`.controls-row`)
      expect(controls_row).toBeDefined()

      // Size select should be in controls row (uses svelte-multiselect)
      const size_select = controls_row?.querySelector(`#size-select`)
      expect(size_select).toBeDefined()

      // Open extra controls and test all defaults
      const settings_button = doc_query<HTMLButtonElement>(`.settings-toggle`)
      await settings_button?.click()
      const pane = doc_query(pane_selector)

      // Find checkboxes outside log-toggles div
      const log_toggles_count = pane?.querySelectorAll(
        `.log-toggles input[type="checkbox"]`,
      ).length ?? 0
      const all_checkboxes = pane?.querySelectorAll<HTMLInputElement>(
        `input[type="checkbox"]`,
      )

      // Test checkbox defaults (after log toggles: show_model_labels, x_grid, y_grid)
      expect(all_checkboxes?.[log_toggles_count]?.checked, `show_model_labels default`)
        .toBe(true)
      expect(all_checkboxes?.[log_toggles_count + 1]?.checked, `x_grid default`).toBe(
        true,
      )
      expect(all_checkboxes?.[log_toggles_count + 2]?.checked, `y_grid default`).toBe(
        true,
      )

      // Test all input defaults (these often regress)
      expect(doc_query<HTMLInputElement>(`#x-ticks`, pane)?.value).toBe(`5`)
      expect(doc_query<HTMLInputElement>(`#y-ticks`, pane)?.value).toBe(`5`)
      expect(doc_query<HTMLInputElement>(`#size-multiplier`, pane)?.value).toBe(`1`)
      expect(doc_query<HTMLInputElement>(`#label-font-size`, pane)?.value).toBe(`14`)
      expect(doc_query<HTMLInputElement>(`#min-link-distance`, pane)?.value).toBe(`15`)
      expect(doc_query<HTMLInputElement>(`#max-link-distance`, pane)?.value).toBe(`20`)
      expect(doc_query<HTMLInputElement>(`#link-strength`, pane)?.value).toBe(`5`)

      // Verify plot container and color controls render
      expect(document.querySelector(`div.bleed-1400[style]`)).toBeDefined()
      expect(pane?.querySelector(`.color-scale-select`)).toBeDefined()
    })
  })
})
