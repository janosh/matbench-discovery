import kappa_103_analysis from '$figs/kappa-103-analysis.jsonl'
import { MODELS } from '$lib'
import { ALL_METRICS } from '$lib/labels'
import { assemble_row_data, get_nested_number } from '$lib/metrics'
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

// mirrors the page's leaderboard filter: phonon metrics + discovery data (energy-only
// models are hidden by MetricsTable's defaults)
const phonon_leaderboard_count = MODELS.filter(
  (model) =>
    get_nested_number(model, `metrics.phonons.kappa_103.κ_SRME`) != null &&
    get_nested_number(model, `metrics.phonons.kappa_103.κ_SRE`) != null &&
    model.targets !== `E` &&
    typeof model.metrics?.discovery === `object` &&
    model.metrics.discovery.unique_prototypes,
).length

const get_headers = (root: ParentNode) =>
  Array.from(root.querySelectorAll(`th`), (header) =>
    header.textContent?.replace(/\s*[↑↓]\s*$/, ``).trim(),
  )

const kappa_sort_select = (): HTMLSelectElement =>
  doc_query<HTMLSelectElement>(`.kappa-model-select select`)

// the model picker is a svelte-multiselect; its selected chip lives in this ul
const kappa_selected_model = (): string | null | undefined =>
  document.querySelector(`.kappa-model-select ul[aria-label="selected options"]`)
    ?.textContent

// open the multiselect dropdown and return its option labels in listed order
async function kappa_model_option_labels(): Promise<string[]> {
  doc_query<HTMLInputElement>(`.kappa-model-select input`).click()
  await tick()
  return [
    ...document.querySelectorAll<HTMLElement>(
      `.kappa-model-select ul[role="listbox"] > li`,
    ),
  ].map((option) => option.textContent?.trim() ?? ``)
}

const heading_texts = (): (string | undefined)[] =>
  [...document.querySelectorAll(`h2`)].map((heading) =>
    heading.textContent?.replaceAll(/\s+/g, ` `).trim(),
  )

describe(`Phonons Task Page`, () => {
  it(`renders page structure, scatter, and diagnostics`, () => {
    mount(PhononsPage, { target: document.body })

    expect(document.querySelector(`h1`)?.textContent).toContain(
      `MLFF Phonon Modeling Metrics`,
    )

    const leaderboard = doc_query(`section.full-bleed:not(.robustness-table)`)
    expect(leaderboard.querySelectorAll(`tbody tr`)).toHaveLength(
      phonon_leaderboard_count,
    )

    const headings = heading_texts()
    expect(headings).toContain(`Failure Modes`)
    expect(headings).toContain(`κSRME vs Params`)

    const scatter = doc_query<HTMLDivElement>(`div.scatter`)
    expect(scatter.getAttribute(`style`)).toContain(`height: 800px`)

    // SRME-vs-kappa scatter and frequency parity plot side by side.
    expect(doc_query(`.diagnostics-grid`).querySelectorAll(`div.scatter`)).toHaveLength(2)
  })

  it(`renders the failure-mode robustness table with one row per kappa model`, () => {
    mount(PhononsPage, { target: document.body })

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
    expect(kappa_col_indices.every((col_idx) => col_idx >= 0)).toBe(true)
    const rows = [...leaderboard.querySelectorAll<HTMLTableRowElement>(`tbody tr`)]
    expect(rows).toHaveLength(phonon_leaderboard_count)
    for (const row of rows) {
      const cells = [...row.querySelectorAll(`td`)]
      for (const col_idx of kappa_col_indices) {
        expect(cells[col_idx].textContent?.trim()).not.toBe(`n/a`)
      }
    }
  })

  it(`κ_SRE column carries κ_SRE values (not κ_SRME)`, () => {
    const rows = assemble_row_data(
      `unique_prototypes`,
      () => true, // model_filter: keep all
      false, // show_energy_only
    )
    const mace = rows.find((row) => String(row.Model).includes(`MACE-MP-0`))
    // mace-mp-0.yml has κ_SRE 0.471 (and κ_SRME 0.6823) — the column must use κ_SRE
    expect(mace?.[ALL_METRICS.κ_SRE.key]).toBe(0.471)
  })

  it(`offers sort controls and defaults the Compare dropdown to κSRME order`, async () => {
    mount(PhononsPage, { target: document.body })

    // sort control offers alphabetical, κSRME, and date-added modes
    expect([...kappa_sort_select().options].map((opt) => opt.value)).toEqual([
      `kappa`,
      `name`,
      `date`,
    ])

    // default order is κSRME ascending (best phonon models first); option labels are
    // `<model_name> (<κSRME>)`, so strip the suffix to look up each model's metric
    const srme_path = `metrics.phonons.kappa_103.κ_SRME`
    const labels = await kappa_model_option_labels()
    const srmes = labels.map((label) => {
      const model_name = label.replace(/ \([^()]*\)$/, ``)
      const model = MODELS.find((mdl) => mdl.model_name === model_name)
      return model ? (get_nested_number(model, srme_path) ?? Infinity) : Infinity
    })
    expect(srmes.filter(Number.isFinite).length).toBeGreaterThan(1)
    expect(srmes).toEqual([...srmes].toSorted((s1, s2) => s1 - s2))
  })

  it(`restores URL state for diagnostics controls`, async () => {
    await mount_with_url(
      PhononsPage,
      `http://localhost/tasks/phonons?model=mace-mp-0&model_sort=name&x=F1&y=κ_SRE`,
    )

    expect(kappa_sort_select().value).toBe(`name`)
    const mace_name = MODELS.find((mdl) => mdl.model_key === `mace-mp-0`)?.model_name
    if (!mace_name) throw new Error(`mace-mp-0 not found in MODELS`)
    expect(kappa_selected_model()).toContain(mace_name)
    expect(heading_texts()).toContain(`κSRE vs F1`)
  })

  it(`restores URL state for table sort and filters`, async () => {
    await mount_with_url(
      PhononsPage,
      `http://localhost/tasks/phonons?sort=κ_SRE&dir=asc&openness=OSOD&heatmap=0`,
    )

    const header = sorted_header()
    expect(header?.textContent).toContain(`κSRE`)
    expect(header?.getAttribute(`aria-sort`)).toBe(`ascending`)
    expect(filter_summary_badge(`Openness`)).toContain(`(1/4)`)
    expect(checkbox_for(`Heatmap`).checked).toBe(false)
  })
})
