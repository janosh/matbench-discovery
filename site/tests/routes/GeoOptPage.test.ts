import spg_sankeys from '$figs/spg-sankeys.jsonl'
import struct_rmsd_cdf from '$figs/struct-rmsd-cdf.jsonl'
import sym_ops_diff from '$figs/sym-ops-diff-bar.jsonl'
import { by_benchmark_added_desc, MODELS } from '$lib'
import GeoOptPage from '$routes/tasks/geo-opt/+page.svelte'
import { describe, expect, it } from 'vitest'
import {
  checkbox_for,
  doc_query,
  filter_summary_badge,
  mount_with_url,
  sorted_header,
} from '../index'

const key_label_pairs = MODELS.flatMap((model) =>
  model.model_key ? [[model.model_key, model.model_name] as const] : [],
)
const label_by_model_key = new Map(key_label_pairs)

const selected_text = (): string | undefined =>
  document.querySelector(`.plot-controls .multiselect ul[aria-label="selected options"]`)
    ?.textContent

const cdf_labels = (): string[] =>
  (document.querySelector(`.rmsd-cdf`)?.getAttribute(`aria-label`) ?? ``)
    .replace(`RMSD CDF models: `, ``)
    .split(`, `)
    .filter(Boolean)

const histogram_labels = (): string[] =>
  [...document.querySelectorAll(`.sym-ops-list figcaption`)].map(
    (caption) => caption.textContent?.replace(/\s+\(σ=.*$/, ``) ?? ``,
  )

const sankey_labels = (): string[] =>
  [...document.querySelectorAll(`.spg-sankeys h3`)].map(
    (heading) => heading.textContent ?? ``,
  )

describe(`Geo Opt Task Page`, () => {
  it(`renders intro, leaderboard, comparison, and diagnostics in order`, async () => {
    await mount_with_url(GeoOptPage, `http://localhost/tasks/geo-opt`)

    expect(doc_query(`h1`).textContent).toContain(`MLFF Geometry Optimization`)
    const section_headings = [...document.querySelectorAll(`h2`)].map((heading) =>
      heading.textContent?.trim(),
    )
    expect(section_headings).toStrictEqual([
      `Leaderboard`,
      `Model Comparison`,
      `Aggregate Diagnostics`,
    ])
    expect(document.body.textContent).toContain(`RMSD is symprec-invariant`)
    expect(doc_query(`.collapsible-legend .scatter`)).toBeInstanceOf(HTMLElement)
    expect(cdf_labels().length).toBeGreaterThan(0)
  })

  // default selection = the 5 most recently added models among those with plot payloads
  it(`preselects the newest models by benchmark_added`, async () => {
    await mount_with_url(GeoOptPage, `http://localhost/tasks/geo-opt`)

    const payload_keys = new Set(
      [...struct_rmsd_cdf.models, ...sym_ops_diff.models, ...spg_sankeys.models].map(
        (model) => model.model_key,
      ),
    )
    // compare dates, not names: models added on the same day tie and may swap order
    const expected_dates = MODELS.filter((model) => payload_keys.has(model.model_key))
      .toSorted(by_benchmark_added_desc)
      .slice(0, 5)
      .map((model) => model.dates.benchmark_added)
    const selected_dates = [
      ...document.querySelectorAll(
        `.plot-controls .multiselect ul[aria-label="selected options"] > li`,
      ),
    ].map((item) => {
      const name = item.textContent?.trim()
      const model = MODELS.find((candidate) => candidate.model_name === name)
      if (!model) throw new Error(`unknown selected model ${name}`)
      return model.dates.benchmark_added
    })
    expect(selected_dates).toStrictEqual(expected_dates)
    expect(new Set(selected_dates).size).toBeGreaterThan(1)
  })

  it(`filters every aggregate plot from the models query param`, async () => {
    const shared_model = spg_sankeys.models.find(
      (model) =>
        struct_rmsd_cdf.models.some((entry) => entry.model_key === model.model_key) &&
        sym_ops_diff.models.some((entry) => entry.model_key === model.model_key),
    )
    if (!shared_model) throw new Error(`No model shared by all geo-opt payloads`)
    const { model_key, label } = shared_model
    await mount_with_url(GeoOptPage, `http://localhost/tasks/geo-opt?models=${model_key}`)

    expect(selected_text()).toContain(label_by_model_key.get(model_key) ?? label)
    expect(cdf_labels()).toStrictEqual([label])
    expect(histogram_labels()).toStrictEqual([label])
    expect(sankey_labels()).toStrictEqual([label])
  })

  it(
    `keeps empty states for unselected aggregate diagnostics`,
    { timeout: 30_000 },
    async () => {
      await mount_with_url(GeoOptPage, `http://localhost/tasks/geo-opt?models=`)

      expect(document.querySelectorAll(`.empty-note`)).toHaveLength(3)
      expect(document.querySelector(`.rmsd-cdf`)).toBeNull()
    },
  )

  it(`restores scatter, sort, and metrics-table filters from URL params`, async () => {
    await mount_with_url(
      GeoOptPage,
      `http://localhost/tasks/geo-opt?x=model_params&y=symmetry_match_1e-5&sort=Model&dir=desc&train=MPtrj&openness=OSOD,OSCD&heatmap=0`,
    )

    const scatter_heading = [...document.querySelectorAll(`h3`)].find((heading) =>
      heading.textContent?.includes(` vs `),
    )
    expect(scatter_heading?.textContent).toContain(`Params`)
    expect(scatter_heading?.textContent).toContain(`Σ`)
    expect(sorted_header()?.textContent).toContain(`Model`)
    expect(sorted_header()?.getAttribute(`aria-sort`)).toBe(`descending`)
    expect(filter_summary_badge(`Training data`)).toContain(`(1)`)
    expect(filter_summary_badge(`Openness`)).toContain(`(2/4)`)
    expect(checkbox_for(`Heatmap`).checked).toBe(false)
  })
})
