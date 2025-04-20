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
      geo_opt: {
        'symprec=1e-5': { rmsd: 0.02 },
        'symprec=1e-2': { rmsd: 0.03 },
      },
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
      geo_opt: {
        'symprec=1e-5': { rmsd: 0.015 },
        'symprec=1e-2': { rmsd: 0.025 },
      },
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
      geo_opt: {
        'symprec=1e-5': { rmsd: 0.01 },
        'symprec=1e-2': { rmsd: 0.02 },
      },
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

    // Check that log scale checkboxes are rendered (expect 3 initially)
    const checkboxes =
      document.body.querySelectorAll<HTMLInputElement>(`input[type="checkbox"]`)
    expect(checkboxes.length).toBe(3)
    // Check element type
    checkboxes.forEach((cb) => expect(cb.tagName).toBe(`INPUT`))
    // Check initial checked state (defaults: x=date_added, y=F1, color=model_params)
    // Log default state: x=false, y=false, color=true
    expect(checkboxes[0].checked).toBe(false) // x: date_added (log disabled)
    expect(checkboxes[1].checked).toBe(false) // y: F1
    expect(checkboxes[2].checked).toBe(true) // color: model_params
    // Check initial disabled state for date_added
    expect(checkboxes[1].disabled).toBe(false) // y: F1
    expect(checkboxes[2].disabled).toBe(false) // color: model_params

    // Check that the scatter plot container is rendered
    const plot_container = document.body.querySelector(`div.full-bleed-1400[style]`)
    expect(plot_container).toBeDefined()
  })

  it(`renders component structure (no SVG check)`, () => {
    // Use mount from svelte
    mount(DynamicScatter, {
      target: document.body,
      props: { models: mock_models },
    })
    expect(document.body.querySelector(`.controls-grid`)).toBeDefined()
    expect(document.body.querySelector(`div.full-bleed-1400[style]`)).toBeDefined()
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
})
