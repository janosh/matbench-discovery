import elem_prev from '$figs/element-prevalence-vs-error.jsonl'
import fp_diff from '$figs/scatter-largest-fp-diff-each-error.jsonl'
import DiscoveryTmiPage from '$routes/tasks/discovery/tmi/+page.svelte'
import { describe, expect, it } from 'vitest'
import { mount_with_url } from '../index'

// the page's per-figure model selects are wrapped in <label>s (the periodic-table
// heatmap's model select is not, which keeps it out of this list)
const selected_texts = (): string[] =>
  [
    ...document.querySelectorAll(`label .multiselect ul[aria-label="selected options"]`),
  ].map((list) => list.textContent ?? ``)

describe(`Discovery TMI Page`, () => {
  it.each([
    [`URL has no model params`, ``],
    [`URL has unknown model tokens`, `?models=bogus&fp_model=bogus`],
  ])(`defaults model selections when %s`, async (_case_name, query) => {
    await mount_with_url(DiscoveryTmiPage, `http://localhost/tasks/discovery/tmi${query}`)

    const [elem_prev_text, fp_text] = selected_texts()
    for (const model of elem_prev.models.slice(0, 3)) {
      expect(elem_prev_text).toContain(model.label)
    }
    expect(fp_text).toContain(fp_diff.models[0].label)
  })

  it(`restores model selections from URL params`, async () => {
    const elem_model = elem_prev.models.at(-1)
    const fp_model = fp_diff.models.at(-1)
    if (!elem_model || !fp_model) throw new Error(`missing payload models`)
    const query = new URLSearchParams({
      models: elem_model.model_key,
      fp_model: fp_model.model_key,
    })
    const url = `http://localhost/tasks/discovery/tmi?${query}`
    await mount_with_url(DiscoveryTmiPage, url)

    const [elem_prev_text, fp_text] = selected_texts()
    expect(elem_prev_text).toContain(elem_model.label)
    // multi-select restored to exactly one model, not the 3 defaults
    for (const model of elem_prev.models.slice(0, 3)) {
      if (model.model_key !== elem_model.model_key) {
        expect(elem_prev_text).not.toContain(model.label)
      }
    }
    expect(fp_text).toContain(fp_model.label)
  })
})
