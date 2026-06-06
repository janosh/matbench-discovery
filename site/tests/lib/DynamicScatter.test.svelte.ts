import DynamicScatter from '$lib/plot/DynamicScatter.svelte'
import type { ModelData } from '$lib/types'
import { mount, tick } from 'svelte'
import { describe, expect, it } from 'vitest'

// Models carry values for the default axes (x=κ_SRME, y=CPS, color=F1, size=Params)
const mock_models: ModelData[] = [
  {
    model_name: `Test Model 1`,
    model_key: `test-model-1`,
    model_params: 100_000,
    date_added: `2023-01-01`,
    color: `blue`,
    CPS: 0.6,
    metrics: {
      discovery: { unique_prototypes: { F1: 0.85 } },
      phonons: { kappa_103: { κ_SRME: 0.4 } },
    },
  } as unknown as ModelData,
  {
    model_name: `Test Model 2`,
    model_key: `test-model-2`,
    model_params: 500_000,
    date_added: `2023-03-15`,
    color: `red`,
    CPS: 0.7,
    metrics: {
      discovery: { unique_prototypes: { F1: 0.9 } },
      phonons: { kappa_103: { κ_SRME: 0.3 } },
    },
  } as unknown as ModelData,
]

// Each model renders as its own matterviz series, so the built-in legend handles per-model
// colors, marker shapes, and toggling. That output lives inside the SVG and only renders
// with real layout, so the legend/markers/colors are verified in the browser, not here.
describe(`DynamicScatter.svelte`, () => {
  it(`mounts with the controls row and a size select defaulting to Params`, async () => {
    // pass a model_filter to exercise that path; its effect on the SVG isn't observable here
    mount(DynamicScatter, {
      target: document.body,
      props: {
        models: mock_models,
        model_filter: (model: ModelData) => model.model_name !== `Test Model 2`,
      },
    })
    await tick()

    const size_select = document.querySelector(
      `div.bleed-1400 .controls-row #size-select`,
    )
    expect(size_select).not.toBeNull()
    expect(size_select?.closest(`.multiselect`)?.textContent).toContain(`Params`)
  })
})
