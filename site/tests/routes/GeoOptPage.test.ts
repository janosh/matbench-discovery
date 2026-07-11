import spg_sankeys from '$figs/spg-sankeys.jsonl'
import struct_rmsd_cdf from '$figs/struct-rmsd-cdf.jsonl'
import sym_ops_diff from '$figs/sym-ops-diff-bar.jsonl'
import { MODELS } from '$lib'
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

  it(`filters every aggregate plot from the models query param`, async () => {
    const shared_model = spg_sankeys.models.find(
      (model) =>
        struct_rmsd_cdf.models.some((entry) => entry.label === model.label) &&
        sym_ops_diff.models.some((entry) => entry.label === model.label),
    )
    if (!shared_model) throw new Error(`No model shared by all geo-opt payloads`)
    const { key, label } = shared_model
    await mount_with_url(GeoOptPage, `http://localhost/tasks/geo-opt?models=${key}`)

    expect(selected_text()).toContain(label_by_model_key.get(key) ?? label)
    expect(cdf_labels()).toStrictEqual([label])
    expect(histogram_labels()).toStrictEqual([label])
    expect(sankey_labels()).toStrictEqual([label])
  })

  it(`keeps empty states for unselected aggregate diagnostics`, async () => {
    await mount_with_url(GeoOptPage, `http://localhost/tasks/geo-opt?models=`)

    expect(document.querySelectorAll(`.empty-note`)).toHaveLength(3)
    expect(document.querySelector(`.rmsd-cdf`)).toBeNull()
  })

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
