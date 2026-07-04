import spg_sankeys from '$figs/spg-sankeys.jsonl'
import sym_ops_diff from '$figs/sym-ops-diff-bar.jsonl'
import { MODELS } from '$lib'
import GeoOptPage from '$routes/tasks/geo-opt/+page.svelte'
import { describe, expect, it } from 'vitest'
import { mount_with_url } from '../index'

const key_label_pairs = MODELS.flatMap((model) =>
  model.model_key ? [[model.model_key, model.model_name] as const] : [],
)
const label_by_model_key = new Map(key_label_pairs)
const model_key_by_label = new Map(key_label_pairs.map(([key, label]) => [label, key]))
const date_added_by_key = new Map(
  MODELS.map((model) => [model.model_key, Date.parse(model.date_added ?? ``) || 0]),
)
// mirror the page's default derivation exactly (same initial order + stable sort) so
// date_added ties at the rank-5 cutoff resolve identically here and in the page
const selectable_model_keys = [
  ...new Set([
    ...sym_ops_diff.models.map(({ label }) => model_key_by_label.get(label) ?? label),
    ...spg_sankeys.models.map(({ key }) => key),
  ]),
]
const default_model_keys = selectable_model_keys
  .toSorted(
    (key_1, key_2) =>
      (date_added_by_key.get(key_2) ?? 0) - (date_added_by_key.get(key_1) ?? 0),
  )
  .slice(0, 5)
const last_default_key = default_model_keys.at(-1) ?? ``

const sort_labels = (labels: string[]): string[] =>
  labels.toSorted((label_1, label_2) => label_1.localeCompare(label_2))

const selected_text = (): string | undefined =>
  document.querySelector(`.multiselect ul[aria-label="selected options"]`)?.textContent

const histogram_labels = (): string[] =>
  [...document.querySelectorAll(`.sym-ops-list figcaption`)].map(
    (caption) => caption.textContent?.replace(/\s+\(σ=.*$/, ``) ?? ``,
  )

const sankey_labels = (): string[] =>
  [...document.querySelectorAll(`.spg-sankeys h3`)].map(
    (heading) => heading.textContent ?? ``,
  )

describe(`Geo Opt Task Page`, () => {
  it.each([
    {
      name: `defaults symmetry plots to the 5 most recently added models`,
      query: ``,
      expected_keys: default_model_keys,
    },
    {
      name: `filters symmetry histograms and sankeys from the models query param`,
      query: `?models=${last_default_key}`,
      expected_keys: [last_default_key],
    },
  ])(`$name`, async ({ query, expected_keys }) => {
    await mount_with_url(GeoOptPage, `http://localhost/tasks/geo-opt${query}`)

    const expected_labels = expected_keys.map(
      (model_key) => label_by_model_key.get(model_key) ?? model_key,
    )
    expect(sort_labels(histogram_labels())).toStrictEqual(sort_labels(expected_labels))
    expect(sort_labels(sankey_labels())).toStrictEqual(sort_labels(expected_labels))
    for (const label of expected_labels) expect(selected_text()).toContain(label)
  })
})
