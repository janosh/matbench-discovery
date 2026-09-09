import elem_prev from '$figs/element-prevalence-vs-error.jsonl'
import fp_diff from '$figs/scatter-largest-fp-diff-each-error.jsonl'
import each_errors from '$figs/scatter-largest-each-errors-fp-diff.jsonl'
import hist_largest from '$figs/hist-largest-each-errors-fp-diff.jsonl'
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
  it.each([null, ``, `bogus`, `constructor`])(
    `defaults unknown model selections: %s`,
    async (value) => {
      const query =
        value === null
          ? ``
          : `?models=${value}&fp_model=${value}&each_model=${value}&hist_model=${value}`
      await mount_with_url(
        DiscoveryTmiPage,
        `http://localhost/tasks/discovery/tmi${query}`,
      )

      const [elem_prev_text, ...single_texts] = selected_texts()
      for (const model of elem_prev.models.slice(0, 3)) {
        expect(elem_prev_text).toContain(model.label)
      }
      for (const [idx, { models }] of [fp_diff, each_errors, hist_largest].entries()) {
        expect(single_texts[idx]).toContain(models[0].label)
      }
    },
  )

  it(`restores model selections from URL params`, async () => {
    const elem_model = elem_prev.models.at(-1)
    const single_models = [fp_diff, each_errors, hist_largest].map(
      ({ models }) => models[models.length - 1],
    )
    if (!elem_model) throw new Error(`missing element prevalence models`)
    const query = new URLSearchParams({
      models: `unknown,${elem_model.model_key},${elem_model.model_key}`,
      fp_model: single_models[0].model_key,
      each_model: single_models[1].model_key,
      hist_model: single_models[2].model_key,
    })
    const url = `http://localhost/tasks/discovery/tmi?${query}`
    await mount_with_url(DiscoveryTmiPage, url)

    const [elem_prev_text, ...single_texts] = selected_texts()
    expect(elem_prev_text).toContain(elem_model.label)
    // multi-select restored to exactly one model, not the 3 defaults
    for (const model of elem_prev.models.slice(0, 3)) {
      if (model.model_key !== elem_model.model_key) {
        expect(elem_prev_text).not.toContain(model.label)
      }
    }
    for (const [idx, model] of single_models.entries()) {
      expect(single_texts[idx]).toContain(model.label)
    }
    expect(new URL(location.href).searchParams.get(`models`)).toBe(elem_model.model_key)
  })
})
