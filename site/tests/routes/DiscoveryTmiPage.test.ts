import elem_prev from '#figs/element-prevalence-vs-error.jsonl'
import hist_clf from '#figs/hist-clf-pred-hull-dist.jsonl'
import DiscoveryTmiPage from '#routes/benchmarks/discovery/tmi/+page.svelte'
import { tick } from 'svelte'
import { describe, expect, it } from 'vite-plus/test'
import { doc_query, mount_with_url } from '../index.js'

const ranked_histograms = hist_clf.models.toSorted((left, right) => right.f1 - left.f1)
const histogram_titles = () =>
  [...document.querySelectorAll(`.histogram-panels figcaption`)].map(
    (caption) => caption.textContent,
  )
const histogram_title = (model: (typeof hist_clf.models)[number]) =>
  `${model.label} · F1=${model.f1}`

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
      const query = value === null ? `` : `?models=${value}&hist_models=${value}`
      await mount_with_url(
        DiscoveryTmiPage,
        `http://localhost/benchmarks/discovery/tmi${query}`,
      )

      const [elem_prev_text] = selected_texts()
      expect(selected_texts()).toHaveLength(1)
      for (const model of elem_prev.models.slice(0, 3)) {
        expect(elem_prev_text).toContain(model.label)
      }
      expect(histogram_titles()).toEqual(
        ranked_histograms.slice(0, 4).map(histogram_title),
      )
    },
  )

  it(`restores model selections from URL params`, async () => {
    const elem_model = elem_prev.models.at(-1)
    if (!elem_model) throw new Error(`missing element prevalence models`)
    const query = new URLSearchParams({
      models: `unknown,${elem_model.model_key},${elem_model.model_key}`,
      hist_models: [...ranked_histograms.slice(-2).toReversed(), ranked_histograms.at(-1)]
        .map((model) => model?.model_key)
        .join(`,`),
    })
    const url = `http://localhost/benchmarks/discovery/tmi?${query}`
    await mount_with_url(DiscoveryTmiPage, url)

    const [elem_prev_text] = selected_texts()
    expect(selected_texts()).toHaveLength(1)
    expect(elem_prev_text).toContain(elem_model.label)
    // multi-select restored to exactly one model, not the 3 defaults
    for (const model of elem_prev.models.slice(0, 3)) {
      if (model.model_key !== elem_model.model_key) {
        expect(elem_prev_text).not.toContain(model.label)
      }
    }
    expect(new URL(location.href).searchParams.get(`models`)).toBe(elem_model.model_key)
    expect(histogram_titles()).toEqual(ranked_histograms.slice(-2).map(histogram_title))
    doc_query<HTMLButtonElement>(`.histogram-selector ul.selected button.remove`).click()
    await tick()
    expect(histogram_titles()).toEqual([
      histogram_title(ranked_histograms[ranked_histograms.length - 1]),
    ])
    expect(new URL(location.href).searchParams.get(`hist_models`)).toBe(
      ranked_histograms.at(-1)?.model_key,
    )
    expect(new URL(location.href).searchParams.get(`models`)).toBe(elem_model.model_key)
  })

  it(`restores more than four histograms, ordered by F1`, async () => {
    const models = ranked_histograms.slice(2, 8)
    await mount_with_url(
      DiscoveryTmiPage,
      `http://localhost/benchmarks/discovery/tmi?hist_models=${models
        .toReversed()
        .map(({ model_key }) => model_key)
        .join(`,`)}`,
    )
    expect(histogram_titles()).toEqual(models.map(histogram_title))
  })
})
