import { MODELS } from '$lib'
import ElementErrorsPtableHeatmap from '$routes/tasks/discovery/tmi/ElementErrorsPtableHeatmap.svelte'
import { per_element_each_errors as per_elem_each_errors } from '$lib/per-element-errors'
import { describe, expect, it } from 'vitest'
import { mount } from '../index'

describe(`ElementErrorsPtableHeatmap`, () => {
  it(`defaults to a model with per-element error data`, () => {
    mount(ElementErrorsPtableHeatmap, { target: document.body })

    // the ModelSelect's single selected chip names the default model (its remove button
    // is icon-only, so the chip text is exactly the model name)
    const chips = document.querySelectorAll(`ul[aria-label="selected options"] li`)
    expect(chips).toHaveLength(1)
    const selected_name = chips[0].textContent?.trim()
    const model = MODELS.find((candidate) => candidate.model_name === selected_name)
    if (!model) throw new Error(`no model named like selected chip ${selected_name}`)

    expect(per_elem_each_errors).toHaveProperty(model.model_key)
  })
})
