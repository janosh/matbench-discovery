import { MODELS } from '$lib'
import { ALL_METRICS } from '$lib/labels'
import { assemble_row_data, get_nested_number } from '$lib/metrics'
import PhononsPage from '$routes/tasks/phonons/+page.svelte'
import { describe, expect, it } from 'vitest'
import { doc_query, mount, mount_with_url } from '../index'

const get_headers = (root: ParentNode) =>
  Array.from(root.querySelectorAll(`th`), (header) =>
    header.textContent?.replace(/\s*[↑↓]\s*$/, ``).trim(),
  )

const kappa_selects = (): HTMLSelectElement[] => [
  ...document.querySelectorAll<HTMLSelectElement>(`.kappa-model-select select`),
]

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

    expect(document.querySelector(`section.full-bleed table`)).not.toBeNull()
    expect(document.querySelectorAll(`tbody tr`).length).toBeGreaterThan(0)

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
    expect(section.querySelectorAll(`td a[href^="/models/"]`).length).toBeGreaterThan(30)
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
    const rows = [...leaderboard.querySelectorAll<HTMLTableRowElement>(`tbody tr`)]
    expect(rows.length).toBeGreaterThan(0)
    for (const row of rows) {
      const cells = [...row.querySelectorAll(`td`)]
      for (const col_idx of kappa_col_indices) {
        expect(col_idx).toBeGreaterThanOrEqual(0)
        expect(cells[col_idx]).toBeDefined()
        expect(cells[col_idx].textContent?.trim()).not.toBe(`n/a`)
      }
    }
  })

  it(`κ_SRE column carries κ_SRE values (not κ_SRME)`, () => {
    const rows = assemble_row_data(
      `unique_prototypes`,
      () => true, // model_filter: keep all
      false, // show_energy_only
      true, // show_non_compliant
      true, // show_compliant
    )
    const mace = rows.find((row) => String(row.Model).includes(`MACE-MP-0`))
    // mace-mp-0.yml has κ_SRE 0.471 (and κ_SRME 0.6823) — the column must use κ_SRE
    expect(mace?.[ALL_METRICS.κ_SRE.key]).toBe(0.471)
  })

  it(`offers sort controls and defaults the Compare dropdown to κSRME order`, () => {
    mount(PhononsPage, { target: document.body })

    const [model_select, sort_select] = kappa_selects()
    // sort control offers alphabetical, κSRME, and date-added modes
    expect([...sort_select.options].map((opt) => opt.value)).toEqual([
      `kappa`,
      `name`,
      `date`,
    ])

    // default order is κSRME ascending (best phonon models first)
    const srme_path = `metrics.phonons.kappa_103.κ_SRME`
    const srmes = [...model_select.options].map((opt) => {
      const model = MODELS.find((mdl) => mdl.model_key === opt.value)
      return model ? (get_nested_number(model, srme_path) ?? Infinity) : Infinity
    })
    expect(srmes.length).toBeGreaterThan(1)
    expect(srmes).toEqual([...srmes].toSorted((s1, s2) => s1 - s2))
  })

  it(`restores URL state for diagnostics controls`, async () => {
    await mount_with_url(
      PhononsPage,
      `http://localhost/tasks/phonons?model=mace-mp-0&model_sort=name&x=F1&y=κ_SRE`,
    )

    expect(kappa_selects().map((select) => select.value)).toEqual([`mace-mp-0`, `name`])
    expect(heading_texts()).toContain(`κSRE vs F1`)
  })
})
