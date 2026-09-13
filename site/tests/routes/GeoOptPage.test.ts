import spg_sankeys from '$figs/spg-sankeys.jsonl'
import struct_rmsd_cdf from '$figs/struct-rmsd-cdf.jsonl'
import sym_ops_diff from '$figs/sym-ops-diff-bar.jsonl'
import { by_benchmark_added_desc } from '$lib'
import { ACTIVE_MODELS, make_table_filters, MODELS } from '$lib/models.svelte'
import { GEO_OPT_SYMMETRY_METRICS, HYPERPARAMS } from '$lib/labels'
import type { ModelData } from '$lib/types'
import { tick } from 'svelte'
import GeoOptPage from '$routes/tasks/geo-opt/+page.svelte'
import { describe, expect, it } from 'vitest'
import {
  checkbox_for,
  doc_query,
  filter_summary_badge,
  mount_with_url,
  sorted_header,
} from '../index'

// Mirrors the page's geo-opt presence check and default table filters.
const geo_opt_row_count = (
  matches: (model: ModelData) => boolean = make_table_filters().matches,
) =>
  ACTIVE_MODELS.filter((model) => model.metrics?.geo_opt != null && matches(model)).length

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
    await mount_with_url(
      GeoOptPage,
      `http://localhost/tasks/geo-opt?models=unknown,${model_key},${model_key}`,
    )

    expect(selected_text()).toContain(label_by_model_key.get(model_key) ?? label)
    expect(cdf_labels()).toStrictEqual([label])
    expect(histogram_labels()).toStrictEqual([label])
    expect(sankey_labels()).toStrictEqual([label])
    expect(new URL(location.href).searchParams.get(`models`)).toBe(model_key)
  })

  it.each([``, `unknown`])(
    `keeps empty states for aggregate diagnostics with models=%s`,
    { timeout: 30_000 },
    async (models) => {
      await mount_with_url(GeoOptPage, `http://localhost/tasks/geo-opt?models=${models}`)

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

  it(`renders table with correct structure, columns, groups, and units`, async () => {
    await mount_with_url(GeoOptPage, `http://localhost/tasks/geo-opt?models=`)

    const table = doc_query(`table`)
    doc_query(`thead`, table)
    doc_query(`tbody`, table)

    const headers = [...document.querySelectorAll(`th`)]
    const header_texts = headers.map((header) => header.textContent?.trim())
    const header_html = headers.map((header) => header.innerHTML)

    expect(header_texts).toContain(`Model`)
    expect(header_texts).not.toContain(`#`)

    // RMSD header omits the unit for a concise default column header
    const rmsd_header = headers.find((header) => header.innerHTML.includes(`RMSD`))
    if (!rmsd_header) throw new Error(`RMSD header not found`)
    expect(rmsd_header.innerHTML).not.toContain(`unitless`)
    expect(rmsd_header.innerHTML).not.toContain(`font-weight: 200`)
    // initial sort is RMSD ascending
    expect(rmsd_header.getAttribute(`aria-sort`)).toBe(`ascending`)

    // f_max header includes unit
    const f_max_header = header_html.find((html) => html.includes(`f<sub>max</sub>`))
    expect(f_max_header).toContain(`(${HYPERPARAMS.max_force.unit})`)

    expect(header_texts).toContain(HYPERPARAMS.ase_optimizer.label)
    expect(header_texts).toContain(HYPERPARAMS.max_steps.label)
    expect(header_texts).toContain(HYPERPARAMS.cell_filter.label)

    // Hidden-by-default columns absent
    expect(header_texts).not.toContain(HYPERPARAMS.n_layers.label)
    expect(
      header_html.some((html) =>
        html.includes(HYPERPARAMS.graph_construction_radius.label),
      ),
    ).toBe(false)

    const symmetry_labels = Object.values(GEO_OPT_SYMMETRY_METRICS).map(
      (metric) => metric.label,
    )
    const found_symmetry = symmetry_labels.filter((label) =>
      header_html.some((html) => html.includes(label)),
    ).length
    expect(found_symmetry).toBe(symmetry_labels.length)

    // Group headers for Symmetry and Hyperparams
    const group_texts = [...document.querySelectorAll(`tr.group-header th`)]
      .map((header) => header.textContent?.trim())
      .filter(Boolean)
    expect(group_texts).toContain(`Symmetry`)
    expect(group_texts).toContain(`Hyperparams`)
  })

  it(`renders geo-opt rows and excludes models without geo-opt metrics`, async () => {
    const model_key = `no-geo-opt-regression`
    ACTIVE_MODELS.push({
      model_key,
      model_name: model_key,
      model_version: `test`,
      targets: `EFS_G`,
      training_sets: [],
      n_training_materials: 1,
      n_training_structures: 1,
      model_params: 1,
      dates: { benchmark_added: `2026-06-30` },
      metrics: {
        discovery: { full_test_set: { F1: 0.1 } },
        diatomics: { energy_mae: 1 },
      },
    } as unknown as ModelData)

    try {
      await mount_with_url(GeoOptPage, `http://localhost/tasks/geo-opt?models=`)

      expect(document.body.textContent).not.toContain(model_key)

      const rows = document.querySelectorAll(`tbody tr`)
      expect(rows).toHaveLength(geo_opt_row_count())

      const model_cells = [...document.querySelectorAll(`td[data-col="Model"]`)]
      expect(model_cells).toHaveLength(rows.length)
      expect(model_cells[0]?.querySelector(`a`)?.getAttribute(`href`)).toMatch(
        /^\/models\//,
      )

      const rmsd_cells = [...document.querySelectorAll(`td`)].filter((td) =>
        td.getAttribute(`data-col`)?.includes(`RMSD`),
      )
      expect(
        rmsd_cells.some((cell) => {
          const text = cell.textContent?.trim()
          return text ? !Number.isNaN(Number(text)) : false
        }),
      ).toBe(true)
    } finally {
      ACTIVE_MODELS.pop()
    }
  })

  it(`toggles heatmap colors`, async () => {
    await mount_with_url(GeoOptPage, `http://localhost/tasks/geo-opt?models=`)

    const checkbox = doc_query<HTMLInputElement>(
      `input[type="checkbox"][aria-label="Toggle heatmap colors"]`,
    )
    expect(checkbox.checked).toBe(true)

    checkbox.click()
    await tick()
    expect(checkbox.checked).toBe(false)

    checkbox.click()
    await tick()
    expect(checkbox.checked).toBe(true)
  })

  it(`opens column visibility panel with checkboxes`, async () => {
    await mount_with_url(GeoOptPage, `http://localhost/tasks/geo-opt?models=`)

    const toggle_btn = doc_query(`.column-toggles summary`)
    toggle_btn.click()
    await tick()

    const column_menu = doc_query(`.sections-container`)
    expect(
      [...column_menu.querySelectorAll(`.section-header`)].map((button) =>
        button.textContent?.trim().replace(/^[▶▼]\s*/, ``),
      ),
    ).toStrictEqual([`Symmetry`, `Hyperparams`])
    expect(column_menu.querySelectorAll(`input[type="checkbox"]`).length).toBeGreaterThan(
      Object.keys(GEO_OPT_SYMMETRY_METRICS).length,
    )
  })

  it.each([`RMSD`, `Model`])(`sorts by %s when header is clicked`, async (col_name) => {
    await mount_with_url(GeoOptPage, `http://localhost/tasks/geo-opt?models=`)

    const headers = [...document.querySelectorAll(`th`)]
    const header = headers.find((candidate) =>
      col_name === `Model`
        ? candidate.textContent?.trim() === `Model`
        : candidate.textContent?.includes(col_name),
    )
    if (!header) throw new Error(`${col_name} header not found`)

    const get_order = () =>
      [...document.querySelectorAll(`td[data-col="Model"]`)].map((cell) =>
        cell.textContent?.trim(),
      )

    const initial = get_order()
    header.click()
    await tick()
    const after_click = get_order()
    header.click()
    await tick()
    const after_second = get_order()

    const changed =
      JSON.stringify(initial) !== JSON.stringify(after_click) ||
      JSON.stringify(after_click) !== JSON.stringify(after_second)
    expect(initial).toHaveLength(geo_opt_row_count())
    expect(changed).toBe(true)
  })

  it.each([
    { training: {}, desc: `no filters` },
    { training: { MPtrj: `require` }, desc: `MPtrj-trained only` },
    { training: { OMat24: `exclude` }, desc: `OMat24 excluded` },
  ] as const)(`filters rows with $desc`, async ({ training }) => {
    const filters = make_table_filters()
    filters.training = { ...training } as typeof filters.training
    const train = Object.entries(training)
      .map(([key, mode]) => (mode === `exclude` ? `-${key}` : key))
      .join(`,`)
    await mount_with_url(
      GeoOptPage,
      `http://localhost/tasks/geo-opt?models=&train=${train}`,
    )

    doc_query(`thead`, doc_query(`table`))

    const rows = document.querySelectorAll(`tbody tr[data-row-idx]`)
    expect(rows).toHaveLength(geo_opt_row_count(filters.matches))
    const first_row = rows[0]
    if (first_row) {
      const href = doc_query<HTMLAnchorElement>(
        `a[href^="/models/"]`,
        first_row,
      ).getAttribute(`href`)
      if (train) doc_query(`button.clear-filters`).click()
      await tick()
      expect(new URL(location.href).searchParams.has(`train`)).toBe(false)
      expect(doc_query(`a[href="${href}"]`, doc_query(`tbody`)).closest(`tr`)).toBe(
        first_row,
      )
    } else expect(doc_query(`tbody`).textContent?.trim()).toBe(`No data`)
  })
})
