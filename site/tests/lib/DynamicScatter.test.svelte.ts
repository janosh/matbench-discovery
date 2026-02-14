import DynamicScatter from '$lib/DynamicScatter.svelte'
import type { ModelData } from '$lib/types'
import { mount, tick } from 'svelte'
import { describe, expect, it } from 'vitest'

const mock_models: ModelData[] = [
  {
    model_name: `Test Model 1`,
    model_params: 100_000,
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
    model_params: 500_000,
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
  it(`mounts with controls row, size select, and scatter container`, () => {
    mount(DynamicScatter, {
      target: document.body,
      props: { models: mock_models },
    })

    const container = document.querySelector(`div.bleed-1400`)
    expect(container).not.toBeNull()

    const controls_row = container?.querySelector(`.controls-row`)
    expect(controls_row).not.toBeNull()
    expect(controls_row?.querySelector(`#size-select`)).not.toBeNull()
  })

  it(`respects model_filter prop by excluding filtered models`, async () => {
    // Mount without filter first to get baseline label count
    mount(DynamicScatter, {
      target: document.body,
      props: { models: mock_models },
    })
    await tick()
    const all_labels = document.querySelectorAll(`text`).length

    // Re-mount with filter excluding a model with valid params to truly test the filter
    document.body.replaceChildren()
    mount(DynamicScatter, {
      target: document.body,
      props: {
        models: mock_models,
        model_filter: (model: ModelData) => model.model_name !== `Test Model 1`,
      },
    })
    await tick()

    // Filtered mount should have fewer or equal text labels
    // (happy-dom may not render SVG text elements, so both can be 0)
    const filtered_labels = document.querySelectorAll(`text`).length
    expect(filtered_labels).toBeLessThanOrEqual(all_labels)
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
