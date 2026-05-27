import { MODELS } from '$lib'
import ElementErrorsPtableHeatmap from '$routes/models/tmi/ElementErrorsPtableHeatmap.svelte'
import per_elem_each_errors from '$routes/models/per-element-each-errors.json'
import { mount } from 'svelte'
import { describe, expect, it } from 'vitest'

describe(`ElementErrorsPtableHeatmap`, () => {
  it(`defaults to a model with per-element error data`, () => {
    const component = mount(ElementErrorsPtableHeatmap, { target: document.body })

    const [current_model] = component.snapshot.capture().current_model
    const model = MODELS.find((candidate) => candidate.model_name === current_model)
    if (!model?.model_key) throw new Error(`missing model_key for ${current_model}`)

    expect(per_elem_each_errors).toHaveProperty(model.model_key)
  })
})
