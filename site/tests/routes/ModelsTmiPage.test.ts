import elem_prev from '$figs/element-prevalence-vs-error.jsonl'
import fp_diff from '$figs/scatter-largest-fp-diff-each-error.jsonl'
import ModelsTmiPage from '$routes/models/tmi/+page.svelte'
import { describe, expect, it } from 'vitest'
import { mount_with_url } from '../index'

// the page's per-figure model selects are wrapped in <label>s (the periodic-table
// heatmap's model select is not, which keeps it out of this list)
const selected_texts = (): (string | undefined)[] =>
  [
    ...document.querySelectorAll(`label .multiselect ul[aria-label="selected options"]`),
  ].map((list) => list.textContent ?? undefined)

describe(`Models TMI Page`, () => {
  it(`defaults model selections when URL has no params`, async () => {
    await mount_with_url(ModelsTmiPage, `http://localhost/models/tmi`)

    // element-prevalence multi-select defaults to the top 3 payload models
    const [elem_prev_text, fp_text] = selected_texts()
    for (const model of elem_prev.models.slice(0, 3)) {
      expect(elem_prev_text).toContain(model.label)
    }
    expect(fp_text).toContain(fp_diff.models[0].label)
  })

  it(`restores model selections from URL params`, async () => {
    const elem_model = elem_prev.models.at(-1)?.label ?? ``
    const fp_model = fp_diff.models.at(-1)?.label ?? ``
    const query = `models=${encodeURIComponent(elem_model)}&fp_model=${encodeURIComponent(
      fp_model,
    )}`
    await mount_with_url(ModelsTmiPage, `http://localhost/models/tmi?${query}`)

    const [elem_prev_text, fp_text] = selected_texts()
    expect(elem_prev_text).toContain(elem_model)
    // multi-select restored to exactly one model, not the 3 defaults
    for (const model of elem_prev.models.slice(0, 3)) {
      if (model.label !== elem_model) expect(elem_prev_text).not.toContain(model.label)
    }
    expect(fp_text).toContain(fp_model)
  })

  it(`falls back to defaults for unknown model tokens`, async () => {
    await mount_with_url(
      ModelsTmiPage,
      `http://localhost/models/tmi?models=bogus&fp_model=bogus`,
    )

    const [elem_prev_text, fp_text] = selected_texts()
    for (const model of elem_prev.models.slice(0, 3)) {
      expect(elem_prev_text).toContain(model.label)
    }
    expect(fp_text).toContain(fp_diff.models[0].label)
  })
})
