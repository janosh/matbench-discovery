import { MODELS } from '$lib'
import ElementErrorsPtableHeatmap from '$routes/tasks/discovery/tmi/ElementErrorsPtableHeatmap.svelte'
import { per_element_each_errors as per_elem_each_errors } from '$lib/per-element-errors'
import { describe, expect, it } from 'vitest'
import { doc_query, mount } from '../index'

describe(`ElementErrorsPtableHeatmap`, () => {
  it(`defaults to a model with per-element error data`, () => {
    mount(ElementErrorsPtableHeatmap, { target: document.body })

    // the ModelSelect's selected chip names the default model
    const selected_name = doc_query(
      `ul[aria-label="selected options"]`,
    ).textContent?.trim()
    const model = MODELS.find((candidate) =>
      selected_name?.includes(candidate.model_name),
    )
    if (!model) throw new Error(`no model matches selected chip ${selected_name}`)

    expect(per_elem_each_errors).toHaveProperty(model.model_key)
  })
})
