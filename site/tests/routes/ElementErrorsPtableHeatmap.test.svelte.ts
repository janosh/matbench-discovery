import { MODELS } from '$lib'
import ElementErrorsPtableHeatmap from '$routes/tasks/discovery/tmi/ElementErrorsPtableHeatmap.svelte'
import { per_element_each_errors as per_elem_each_errors } from '$lib/per-element-errors'
import { describe, expect, it } from 'vitest'
import { mount } from '../index'

describe(`ElementErrorsPtableHeatmap`, () => {
  it(`defaults to a model with per-element error data`, () => {
    const component = mount(ElementErrorsPtableHeatmap, { target: document.body })

    const [current_model] = component.snapshot.capture().current_model
    const model = MODELS.find((candidate) => candidate.model_name === current_model)
    if (!model) throw new Error(`missing model for ${current_model}`)

    expect(per_elem_each_errors).toHaveProperty(model.model_key)
  })
})
