import DynamicScatter from '$lib/DynamicScatter.svelte'
import type { ModelData } from '$lib/types'
import { mount, tick } from 'svelte'
import { describe, expect, it } from 'vitest'

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
    model_params: null,
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
  it(`mounts correctly with default props`, () => {
    mount(DynamicScatter, {
      target: document.body,
      props: { models: mock_models },
    })

    // Controls row with marker size select should render
    const controls_row = document.querySelector(`.controls-row`)
    expect(controls_row).not.toBeNull()

    // Scatter plot container should render
    const plot_container = document.querySelector(`div.bleed-1400[style]`)
    expect(plot_container).not.toBeNull()

    // Size select should be present in controls row
    const size_select = controls_row?.querySelector(`#size-select`)
    expect(size_select).not.toBeNull()
  })

  it(`renders scatter plot container with correct structure`, () => {
    mount(DynamicScatter, { target: document.body, props: { models: mock_models } })

    // ScatterPlot renders inside the bleed container
    const container = document.querySelector(`div.bleed-1400`)
    expect(container).not.toBeNull()
    expect(container?.querySelector(`.controls-row`)).not.toBeNull()
  })

  it(`respects model_filter prop`, () => {
    mount(DynamicScatter, {
      target: document.body,
      props: {
        models: mock_models,
        model_filter: (model: ModelData) => model.model_params !== null,
      },
    })
    // Component should still render with filtered models
    expect(document.querySelector(`.controls-row`)).not.toBeNull()
  })

  it(`renders marker size select with correct default`, async () => {
    mount(DynamicScatter, { target: document.body, props: { models: mock_models } })
    await tick()

    // The size select should show "Params" as the default selected option
    const size_select = document.querySelector(`#size-select`)
    expect(size_select).not.toBeNull()
    // The selected value text should contain "Params"
    const selected_text = size_select?.closest(`.multiselect`)?.textContent
    expect(selected_text).toContain(`Params`)
  })
})
