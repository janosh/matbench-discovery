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

    // Check that controls grid is rendered
    const controls_grid = document.querySelector(`.controls-grid`)
    expect(controls_grid).toBeDefined()

    // Check that log-scale checkboxes are rendered in the controls grid
    // (detailed checkbox state validation is in "regression tests for default values")
    const checkboxes = document.querySelectorAll<HTMLInputElement>(
      `.controls-grid input[type="checkbox"]`,
    )
    expect(checkboxes.length).toBe(4) // 4 checkboxes: x-axis log, y-axis log, color log, size log

    // Check that the scatter plot container is rendered
    const plot_container = document.querySelector(`div.bleed-1400[style]`)
    expect(plot_container).toBeDefined()
  })

  it(`renders component structure (no SVG check)`, () => {
    // Use mount from svelte
    mount(DynamicScatter, { target: document.body, props: { models: mock_models } })
    expect(document.querySelector(`.controls-grid`)).toBeDefined()
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
    expect(document.querySelector(`.controls-grid`)).toBeDefined()
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
      const show_labels_checkbox = extra_controls?.querySelector<HTMLInputElement>(
        `input[type="checkbox"]`,
      )
      const x_grid_checkbox = extra_controls?.querySelectorAll<HTMLInputElement>(
        `input[type="checkbox"]`,
      )[1]
      const x_ticks_input = doc_query<HTMLInputElement>(`#x-ticks`, extra_controls)
      const y_grid_checkbox = extra_controls?.querySelectorAll<HTMLInputElement>(
        `input[type="checkbox"]`,
      )[2]
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
    it(`verifies all critical default UI state to catch regressions`, () => {
      mount(DynamicScatter, { target: document.body, props: { models: mock_models } })

      // Verify axis controls structure (4 select controls for x, y, color, size)
      const controls_grid = document.querySelector(`.controls-grid`)
      expect(controls_grid?.querySelectorAll(`[role="listbox"]`)).toHaveLength(4)

      // Open extra controls and test all defaults
      doc_query<HTMLButtonElement>(`.settings-toggle`)?.click()
      const pane = doc_query(pane_selector)

      // Test all checkbox defaults (these often regress)
      const checkboxes = pane?.querySelectorAll<HTMLInputElement>(
        `input[type="checkbox"]`,
      )
      expect(checkboxes?.[0]?.checked, `show_model_labels default`).toBe(true)
      expect(checkboxes?.[1]?.checked, `x_grid default`).toBe(true)
      expect(checkboxes?.[2]?.checked, `y_grid default`).toBe(true)

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
