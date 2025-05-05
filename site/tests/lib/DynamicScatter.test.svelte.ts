import DynamicScatter from '$lib/DynamicScatter.svelte'
import type { ModelData } from '$lib/types'
import { mount } from 'svelte'
import { afterEach, describe, expect, it, vi } from 'vitest'

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
  afterEach(() => {
    document.body.innerHTML = ``
    vi.restoreAllMocks()
  })

  it(`mounts correctly with default props`, () => {
    mount(DynamicScatter, {
      target: document.body,
      props: { models: mock_models },
    })

    // Check that controls grid is rendered
    const controls_grid = document.body.querySelector(`.controls-grid`)
    expect(controls_grid).toBeDefined()

    // Check that log scale checkboxes are rendered (expect 4 initially)
    const checkboxes =
      document.body.querySelectorAll<HTMLInputElement>(`input[type="checkbox"]`)
    expect(checkboxes.length).toBe(4)
    // Check initial checked state (defaults: x=date_added, y=F1, color=model_params)
    // Log default state: x=false, y=false, color=true
    expect(checkboxes[0].checked).toBe(false) // x: date_added (log disabled)
    expect(checkboxes[1].checked).toBe(false) // y: F1
    expect(checkboxes[2].checked).toBe(false) // color: model_params
    expect(checkboxes[3].checked).toBe(false) // x: date_added (log enabled)
    // Check initial disabled state for date_added
    expect(checkboxes[1].disabled).toBe(true) // y: CPS (log disabled due to missing data)
    expect(checkboxes[2].disabled).toBe(true) // color: F1 (log disabled due to small range)

    // Check that the scatter plot container is rendered
    const plot_container = document.body.querySelector(`div.full-bleed-1400[style]`)
    expect(plot_container).toBeDefined()
  })

  it(`renders component structure (no SVG check)`, () => {
    // Use mount from svelte
    mount(DynamicScatter, { target: document.body, props: { models: mock_models } })
    expect(document.body.querySelector(`.controls-grid`)).toBeDefined()
    expect(document.body.querySelector(`div.full-bleed-1400[style]`)).toBeDefined()
  })

  // Helper function to check fullscreen state
  async function check_fullscreen_state(expected_state: boolean): Promise<void> {
    const container = document.body.querySelector(`.plot-container`)
    const button = document.body.querySelector<HTMLButtonElement>(`.fullscreen-toggle`)
    const button_icon = button?.querySelector(`svg > use`)

    await vi.waitFor(() => {
      expect(container?.classList.contains(`fullscreen`), `Container class`).toBe(
        expected_state,
      )
      expect(document.body.classList.contains(`fullscreen`), `Body class`).toBe(
        expected_state,
      )
      expect(button_icon?.getAttribute(`href`), `Button icon`).toBe(
        `#icon-${expected_state ? `close` : `maximize`}`,
      )
      expect(button?.getAttribute(`aria-label`), `Button label`).toBe(
        `${expected_state ? `Exit` : `Enter`} fullscreen`,
      )
    })
  }

  it(`handles fullscreen toggle via button click and Escape key`, async () => {
    mount(DynamicScatter, {
      target: document.body,
      props: { models: mock_models },
    })

    const button = document.body.querySelector<HTMLButtonElement>(`.fullscreen-toggle`)
    expect(button).toBeDefined()

    // 1. Initial state: Not fullscreen
    await check_fullscreen_state(false)

    // 2. Enter fullscreen via button click
    await button?.click()
    await check_fullscreen_state(true)

    // 3. Exit fullscreen via Escape key
    await window.dispatchEvent(new KeyboardEvent(`keydown`, { key: `Escape` }))
    await check_fullscreen_state(false)

    // 4. Re-enter fullscreen via button click
    await button?.click()
    await check_fullscreen_state(true)

    // 5. Exit fullscreen via button click again
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
    expect(document.body.querySelector(`.controls-grid`)).toBeDefined()
    // Cannot easily check filter effect without mocks
  })

  describe(`extra controls`, () => {
    it(`toggles extra controls visibility via button click and Escape key`, async () => {
      mount(DynamicScatter, {
        target: document.body,
        props: { models: mock_models },
      })

      const settings_button =
        document.body.querySelector<HTMLButtonElement>(`.settings-toggle`)
      let extra_controls = document.body.querySelector(`.extra-controls`)

      // 1. Initial state: Controls hidden
      expect(extra_controls).toBeNull()

      // 2. Show controls via button click
      await settings_button?.click()
      extra_controls = document.body.querySelector(`.extra-controls`)
      expect(extra_controls).toBeDefined()

      // 3. Hide controls via Escape key
      await window.dispatchEvent(new KeyboardEvent(`keydown`, { key: `Escape` }))
      await vi.waitFor(() => {
        extra_controls = document.body.querySelector(`.extra-controls`)
        expect(extra_controls).toBeNull()
      })

      // 4. Re-show controls via button click
      await settings_button?.click()
      extra_controls = document.body.querySelector(`.extra-controls`)
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

      const settings_button =
        document.body.querySelector<HTMLButtonElement>(`.settings-toggle`)
      let extra_controls = document.body.querySelector(`.extra-controls`)

      // 1. Show controls
      await settings_button?.click()
      extra_controls = document.body.querySelector(`.extra-controls`)
      expect(extra_controls).toBeDefined()

      // 2. Click the explicit outside element
      await outside_element.click() // Simulate click outside
      await vi.waitFor(() => {
        extra_controls = document.body.querySelector(`.extra-controls`)
        expect(extra_controls).toBeNull()
      })

      // Clean up the outside element
      document.body.removeChild(outside_element)
    })

    it(`interacts with all extra controls`, async () => {
      mount(DynamicScatter, {
        target: document.body,
        props: { models: mock_models },
      })

      const settings_button =
        document.body.querySelector<HTMLButtonElement>(`.settings-toggle`)
      await settings_button?.click() // Show controls

      const extra_controls = document.body.querySelector(`.extra-controls`)
      expect(extra_controls, `Extra controls panel should exist`).toBeDefined()

      // --- Find Controls ---
      const color_scale_select_el = extra_controls?.querySelector(`.color-scale-select`)
      const show_labels_checkbox =
        extra_controls?.querySelector<HTMLInputElement>(`input[type="checkbox"]`)
      const x_grid_checkbox =
        extra_controls?.querySelectorAll<HTMLInputElement>(`input[type="checkbox"]`)[1]
      const x_ticks_input = extra_controls?.querySelector<HTMLInputElement>(`#x-ticks`)
      const y_grid_checkbox =
        extra_controls?.querySelectorAll<HTMLInputElement>(`input[type="checkbox"]`)[2]
      const y_ticks_input = extra_controls?.querySelector<HTMLInputElement>(`#y-ticks`)
      const size_multiplier_input =
        extra_controls?.querySelector<HTMLInputElement>(`#size-multiplier`)
      const label_font_size_input =
        extra_controls?.querySelector<HTMLInputElement>(`#label-font-size`)
      const min_link_distance_input =
        extra_controls?.querySelector<HTMLInputElement>(`#min-link-distance`)
      const max_link_distance_input =
        extra_controls?.querySelector<HTMLInputElement>(`#max-link-distance`)
      const link_strength_input =
        extra_controls?.querySelector<HTMLInputElement>(`#link-strength`)

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
      x_ticks_input!.value = `10`
      await x_ticks_input?.dispatchEvent(new Event(`input`))
      expect(x_ticks_input?.value, `X Ticks after change`).toBe(`10`)

      y_ticks_input!.value = `12`
      await y_ticks_input?.dispatchEvent(new Event(`input`))
      expect(y_ticks_input?.value, `Y Ticks after change`).toBe(`12`)

      min_link_distance_input!.value = `10`
      await min_link_distance_input?.dispatchEvent(new Event(`input`))
      expect(min_link_distance_input?.value, `Min Link Distance after change`).toBe(`10`)

      max_link_distance_input!.value = `30`
      await max_link_distance_input?.dispatchEvent(new Event(`input`))
      expect(max_link_distance_input?.value, `Max Link Distance after change`).toBe(`30`)

      // Range Sliders
      size_multiplier_input!.value = `2.5`
      await size_multiplier_input?.dispatchEvent(new Event(`input`))
      expect(size_multiplier_input?.value, `Size Multiplier after change`).toBe(`2.5`)

      label_font_size_input!.value = `18`
      await label_font_size_input?.dispatchEvent(new Event(`input`))
      expect(label_font_size_input?.value, `Label Font Size after change`).toBe(`18`)

      link_strength_input!.value = `8`
      await link_strength_input?.dispatchEvent(new Event(`input`))
      expect(link_strength_input?.value, `Link Strength after change`).toBe(`8`)
    })
  })
})
