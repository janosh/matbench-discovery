import kappa_103_analysis from '$figs/kappa-103-analysis.jsonl'
import { MODELS } from '$lib'
import { ALL_METRICS } from '$lib/labels'
import { assemble_row_data, get_nested_number, label_data_path } from '$lib/metrics'
import { make_table_filters } from '$lib/models.svelte'
import PhononsPage from '$routes/tasks/phonons/+page.svelte'
import { tick } from 'svelte'
import { describe, expect, it } from 'vitest'
import {
  checkbox_for,
  doc_query,
  filter_summary_badge,
  mount,
  mount_with_url,
  sorted_header,
} from '../index'

const default_filters = make_table_filters()
const kappa_srme_path = label_data_path(ALL_METRICS.κ_SRME)
const kappa_sre_path = label_data_path(ALL_METRICS.κ_SRE)
// Mirrors the page's model_filter plus MetricsTable's default filters.
const phonon_leaderboard_count = MODELS.filter(
  (model) =>
    get_nested_number(model, kappa_srme_path) != null &&
    get_nested_number(model, kappa_sre_path) != null &&
    default_filters.matches(model),
).length

const get_headers = (root: ParentNode) =>
  Array.from(root.querySelectorAll(`th`), (header) =>
    header.textContent?.replace(/\s*[↑↓]\s*$/, ``).trim(),
  )

const kappa_sort_select = (): HTMLSelectElement =>
  doc_query<HTMLSelectElement>(`.kappa-model-select select`)

// the model picker is a svelte-multiselect; its selected chip lives in this ul
const kappa_selected_model = (): string | null =>
  doc_query(`.kappa-model-select ul[aria-label="selected options"]`).textContent

// Open the multiselect dropdown and return its options in listed order.
async function kappa_model_options(): Promise<HTMLElement[]> {
  doc_query<HTMLInputElement>(`.kappa-model-select input`).click()
  await tick()
  return Array.from(
    document.querySelectorAll<HTMLElement>(`.kappa-model-select ul[role="listbox"] > li`),
  )
}

const heading_texts = (): (string | undefined)[] =>
  [...document.querySelectorAll(`h2`)].map((heading) =>
    heading.textContent?.replaceAll(/\s+/g, ` `).trim(),
  )

describe(`Phonons Task Page`, () => {
  it(`renders the task narrative, scatter, and diagnostics`, () => {
    mount(PhononsPage, { target: document.body })

    expect(doc_query(`h1`).textContent).toContain(`MLFF Phonon Modeling Metrics`)
    expect(document.body.textContent).toContain(`wrong microscopic reasons`)

    const headings = heading_texts()
    expect(headings).toContain(`Failure Diagnostics`)
    expect(headings).toContain(`Model Comparison: κSRME vs κSRE`)

    const scatter = doc_query<HTMLDivElement>(`div.scatter`)
    expect(scatter.getAttribute(`style`)).toContain(`height: 800px`)

    // SRME-vs-kappa scatter and frequency parity plot side by side.
    expect(doc_query(`.diagnostics-grid`).querySelectorAll(`div.scatter`)).toHaveLength(2)

    const section = doc_query(`section.robustness-table`)
    const headers = get_headers(section)
    for (const column of [`κ failed`, `Imag. modes`, `Spectrum W1`]) {
      expect(headers, `missing column ${column}`).toContain(column)
    }
    // one model-page link per entry in the kappa-103 analysis payload
    expect(section.querySelectorAll(`td a[href^="/models/"]`)).toHaveLength(
      kappa_103_analysis.models.length,
    )
  })

  it(`shows only phonon leaderboard columns and rows with phonon metrics`, () => {
    mount(PhononsPage, { target: document.body })

    const leaderboard = doc_query(`section.full-bleed`)
    const headers = get_headers(leaderboard)
    for (const column of [`Model`, `Links`, `κSRME`, `κSRE`]) {
      expect(headers).toContain(column)
    }
    for (const column of [`F1`, `DAF`, `Acc`, `Prec`]) {
      expect(headers).not.toContain(column)
    }

    const kappa_col_indices = [`κSRME`, `κSRE`].map((header) => headers.indexOf(header))
    expect(kappa_col_indices.every((column_idx) => column_idx >= 0)).toBe(true)
    const rows = [...leaderboard.querySelectorAll<HTMLTableRowElement>(`tbody tr`)]
    expect(rows).toHaveLength(phonon_leaderboard_count)
    for (const row of rows) {
      const cells = [...row.querySelectorAll(`td`)]
      for (const column_idx of kappa_col_indices) {
        expect(cells[column_idx].textContent?.trim()).not.toBe(`n/a`)
      }
    }
    expect(sorted_header()?.textContent).toContain(`κSRME`)
    expect(sorted_header()?.getAttribute(`aria-sort`)).toBe(`ascending`)
  })

  it(`κ_SRE column carries κ_SRE values (not κ_SRME)`, () => {
    const rows = assemble_row_data(
      `unique_prototypes`,
      () => true, // model_filter: keep all
      () => true, // filter_matches: no training/openness/targets filters
    )
    const mace_row = rows.find((row) => row.Model.includes(`MACE-MP-0`))
    // mace-mp-0.yml has κ_SRE 0.471 (and κ_SRME 0.6823) — the column must use κ_SRE
    expect(mace_row?.[ALL_METRICS.κ_SRE.key]).toBe(0.471)
  })

  it(`orders Compare models by κSRME and updates the model URL`, async () => {
    await mount_with_url(PhononsPage, `http://localhost/tasks/phonons`)

    // sort control offers alphabetical, κSRME, and date-added modes
    expect([...kappa_sort_select().options].map((option) => option.value)).toEqual([
      `kappa`,
      `name`,
      `date`,
    ])

    // default order is κSRME ascending (best phonon models first); option labels are
    // `<model_name> (<κSRME>)`, so strip the suffix to look up each model's metric
    const current_selection = kappa_selected_model() ?? ``
    const options = await kappa_model_options()
    const srme_values = options.map((option) => {
      const label = option.textContent?.trim() ?? ``
      const model_name = label.replace(/ \([^()]*\)$/, ``)
      const model = MODELS.find((candidate) => candidate.model_name === model_name)
      return model ? (get_nested_number(model, kappa_srme_path) ?? Infinity) : Infinity
    })
    expect(srme_values.filter(Number.isFinite).length).toBeGreaterThan(1)
    expect(srme_values).toEqual(
      [...srme_values].toSorted((value_1, value_2) => value_1 - value_2),
    )

    const target_model = [`chgnet-0.3.0`, `mace-mp-0`]
      .map((model_key) => MODELS.find((model) => model.model_key === model_key))
      .find((model) => model && !current_selection.includes(model.model_name))
    if (!target_model) throw new Error(`No alternate kappa model found`)
    const target_option = options.find((option) =>
      option.textContent?.includes(target_model.model_name),
    )
    if (!target_option)
      throw new Error(`No kappa option found for ${target_model.model_name}`)
    target_option.click()
    await tick()

    expect(kappa_selected_model()).toContain(target_model.model_name)
    expect(new URL(location.href).searchParams.get(`model`)).toBe(target_model.model_key)
  })

  it(`restores URL state for diagnostics, table sort, and filters`, async () => {
    const url = `http://localhost/tasks/phonons?model=mace-mp-0&model_sort=name&x=F1&y=κ_SRE&sort=κ_SRE&dir=asc&openness=OSOD&heatmap=0`
    await mount_with_url(PhononsPage, url)

    expect(kappa_sort_select().value).toBe(`name`)
    const mace_name = MODELS.find((model) => model.model_key === `mace-mp-0`)?.model_name
    if (!mace_name) throw new Error(`mace-mp-0 not found in MODELS`)
    expect(kappa_selected_model()).toContain(mace_name)
    expect(heading_texts()).toContain(`Model Comparison: κSRE vs F1`)

    const header = sorted_header()
    expect(header?.textContent).toContain(`κSRE`)
    expect(header?.getAttribute(`aria-sort`)).toBe(`ascending`)
    expect(filter_summary_badge(`Openness`)).toContain(`(1/4)`)
    expect(checkbox_for(`Heatmap`).checked).toBe(false)
  })
})
