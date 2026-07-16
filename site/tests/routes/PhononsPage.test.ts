import kappa_103_analysis from '$figs/kappa-103-analysis.jsonl'
import { ACTIVE_MODELS } from '$lib'
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
// Mirrors the page's model_filter plus MetricsTable's default filters.
const phonon_leaderboard_count = ACTIVE_MODELS.filter(
  (model) => model.metrics?.phonons?.kappa_103 != null && default_filters.matches(model),
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
    expect(document.body.textContent).toContain(`cancellation can hide errors`)
    expect(document.body.textContent).not.toContain(`working on extending the test set`)
    expect(document.body.textContent).toMatch(/censored\s+to 2 because/)

    const headings = heading_texts()
    expect(headings).toContain(`Model Comparison: κSRME vs κSRE`)

    const scatter = doc_query<HTMLDivElement>(`div.scatter`)
    expect(scatter.getAttribute(`style`)).toContain(`height: 800px`)

    // SRME-vs-kappa scatter and frequency parity plot side by side.
    expect(doc_query(`.diagnostics-grid`).querySelectorAll(`div.scatter`)).toHaveLength(2)
  })

  it(`shows only phonon leaderboard columns and rows with phonon metrics`, () => {
    mount(PhononsPage, { target: document.body })

    const leaderboard = doc_query(`section.full-bleed`)
    const headers = get_headers(leaderboard)
    for (const column of [
      `Model`,
      `Links`,
      `κSRME`,
      `κSRE`,
      `κSRD`,
      `κ failed`,
      `Im(ω)`,
      `W1(ω)`,
    ]) {
      expect(headers).toContain(column)
    }
    for (const column of [`F1`, `DAF`, `Acc`, `Prec`]) {
      expect(headers).not.toContain(column)
    }

    const kappa_col_indices = [`κSRME`, `κSRE`, `κSRD`, `κ failed`, `Im(ω)`].map(
      (header) => headers.indexOf(header),
    )
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

  it(`matches payload and canonical aggregate diagnostics`, () => {
    const n_materials = kappa_103_analysis.material_ids.length
    for (const values of [
      kappa_103_analysis.formulas,
      kappa_103_analysis.spg_nums,
      kappa_103_analysis.kappa_dft,
    ]) {
      expect(values).toHaveLength(n_materials)
    }
    expect(
      kappa_103_analysis.spg_nums.every(
        (spg_num) => Number.isInteger(spg_num) && spg_num >= 1 && spg_num <= 230,
      ),
    ).toBe(true)

    for (const entry of kappa_103_analysis.models) {
      for (const values of [
        entry.kappa_ml,
        entry.srme,
        entry.srme_censored,
        entry.imag_modes,
        entry.broken_sym,
        entry.max_steps,
        entry.freq_w1,
      ]) {
        expect(values).toHaveLength(n_materials)
      }
      const model = ACTIVE_MODELS.find((candidate) => candidate.model_key === entry.key)
      const metrics = model?.metrics?.phonons?.kappa_103
      if (!metrics) throw new Error(`missing phonon metrics for ${entry.key}`)
      const payload_failure_rate =
        entry.srme_censored.filter((value) => value === true).length / n_materials
      expect(payload_failure_rate).toBeCloseTo(metrics.κ_failure_rate, 4)

      const payload_imaginary_rate =
        entry.imag_modes.filter((value) => value === true).length / n_materials
      expect(payload_imaginary_rate).toBeCloseTo(metrics.imaginary_mode_rate, 4)
      if (metrics.spectrum_w1 === null) {
        expect(entry.freq_w1_mean).toBeNull()
      } else {
        expect(entry.freq_w1_mean).toBeCloseTo(metrics.spectrum_w1, 4)
      }
    }
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
      const model = ACTIVE_MODELS.find((candidate) => candidate.model_name === model_name)
      return model?.metrics?.phonons?.kappa_103?.κ_SRME ?? Infinity
    })
    expect(srme_values.filter(Number.isFinite).length).toBeGreaterThan(1)
    expect(srme_values).toEqual(
      [...srme_values].toSorted((value_1, value_2) => value_1 - value_2),
    )

    const target_model = [`chgnet-0.3.0`, `mace-mp-0`]
      .map((model_key) => ACTIVE_MODELS.find((model) => model.model_key === model_key))
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
    const mace_name = ACTIVE_MODELS.find(
      (model) => model.model_key === `mace-mp-0`,
    )?.model_name
    if (!mace_name) throw new Error(`mace-mp-0 not found in ACTIVE_MODELS`)
    expect(kappa_selected_model()).toContain(mace_name)
    expect(heading_texts()).toContain(`Model Comparison: κSRE vs F1`)

    const header = sorted_header()
    expect(header?.textContent).toContain(`κSRE`)
    expect(header?.getAttribute(`aria-sort`)).toBe(`ascending`)
    expect(filter_summary_badge(`Openness`)).toContain(`(1/4)`)
    expect(checkbox_for(`Heatmap`).checked).toBe(false)
  })
})
