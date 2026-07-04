import Page from '$routes/+page.svelte'
import { tick } from 'svelte'
import { beforeEach, describe, expect, it } from 'vitest'
import { doc_query, mount, mount_with_url, sorted_header } from '../index'

const header_text = () =>
  [...document.querySelectorAll(`thead th`)].map((th) => th.textContent).join(` `)
const preset_button = (label: string): HTMLButtonElement => {
  const button = [
    ...document.querySelectorAll<HTMLButtonElement>(`.selection-toggle button`),
  ].find((btn) => btn.textContent?.trim() === label)
  if (!button) throw new Error(`No preset button found for ${label}`)
  return button
}
const table_header = (label: string): HTMLTableCellElement => {
  const header = [...document.querySelectorAll<HTMLTableCellElement>(`thead th`)].find(
    (th) => th.textContent?.replace(/\s*[↑↓]\s*$/, ``).trim() === label,
  )
  if (!header) throw new Error(`no table header labeled ${label}`)
  return header
}

describe(`Landing Page`, () => {
  // happy-dom mounts of the full-column metrics table are slow in CI
  beforeEach(() => {
    mount(Page, { target: document.body })
  }, 30_000)

  // find a SelectToggle bar by one of its button labels (order-independent)
  const toggle_with = (label: string): Element => {
    const toggle = [...document.querySelectorAll(`.selection-toggle`)].find((el) =>
      el.textContent?.includes(label),
    )
    if (!toggle) throw new Error(`no .selection-toggle contains ${label}`)
    return toggle
  }
  const toggle_buttons = (label: string): HTMLButtonElement[] => [
    ...toggle_with(label).querySelectorAll<HTMLButtonElement>(`button`),
  ]
  const button_texts = (label: string) =>
    toggle_buttons(label).map((button) => button.textContent?.trim())
  const select_preset = async (label: string) => {
    preset_button(label).click()
    await tick()
  }

  it(`renders column preset + discovery set toggle buttons`, () => {
    // column-preset bar always shows; the test-set bar only shows in Discovery preset
    expect(button_texts(`Phonons`)).toStrictEqual([
      `Discovery`,
      `Phonons`,
      `Geo Opt`,
      `MD`,
      `Diatomics`,
    ])
    expect(button_texts(`Full Test Set`)).toStrictEqual([
      `Full Test Set`,
      `Unique Prototypes`,
      `10k Most Stable`,
    ])
  })

  // each non-Discovery preset reveals one of its signature columns (the marker) that
  // the default Discovery view hides, while hiding the discovery-only DAF column
  it.each([
    [`Phonons`, `SRE`], // κ_SRE (matches the phonons task page)
    [`Geo Opt`, `Σ`], // symmetry metrics (Σ= / Σ↓ / Σ↑)
    [`MD`, `vDOS`], // vDOS err (RDF is hidden from leaderboards as redundant)
    [`Diatomics`, `E jump`],
  ])(`%s preset shows its metrics + headline cols, hides DAF`, async (preset, marker) => {
    expect(header_text()).toContain(`DAF`) // Discovery is the default preset
    expect(header_text()).not.toContain(marker) // preset-specific col hidden by default

    await select_preset(preset)
    expect(header_text()).toContain(marker) // preset-specific column now visible
    expect(header_text()).not.toContain(`DAF`) // discovery-only column hidden
    for (const headline of [`CPS`, `F1`, `RMSD`]) {
      expect(header_text()).toContain(headline) // headline cols persist across presets
    }
  })

  it(`surfaces a beta warning for the MD metrics`, async () => {
    const n_beta_warnings = () =>
      [...document.querySelectorAll(`blockquote`)].filter((bq) =>
        bq.textContent?.includes(`interpret with caution`),
      ).length
    expect(n_beta_warnings()).toBe(1) // always shown in the page's MD note

    await select_preset(`MD`)
    expect(n_beta_warnings()).toBe(2) // + contextual warning above the MD table
  })

  it(`shows the test-set toggle only in the Discovery preset`, async () => {
    const test_set_shown = () =>
      [...document.querySelectorAll(`.selection-toggle`)].some((toggle) =>
        toggle.textContent?.includes(`Full Test Set`),
      )
    expect(test_set_shown()).toBe(true) // Discovery is the default preset

    await select_preset(`MD`)
    expect(test_set_shown()).toBe(false) // hidden outside Discovery

    await select_preset(`Discovery`)
    expect(test_set_shown()).toBe(true)
  })

  it(`auto-sorts presets until the user manually sorts the table`, async () => {
    expect(header_text()).toContain(`CPS ↑`)

    await select_preset(`MD`)
    expect(header_text()).toContain(`CMDS ↑`)

    table_header(`F1`).click()
    await tick()
    expect(header_text()).toMatch(/F1 [↑↓]/)

    await select_preset(`Geo Opt`)
    expect(header_text()).toMatch(/F1 [↑↓]/)
    expect(header_text()).not.toContain(`RMSD ↓`)
  })

  it(`keeps preset auto-sort after non-header table clicks`, async () => {
    doc_query(`section.full-bleed`).click()
    await tick()

    await select_preset(`MD`)
    expect(header_text()).toContain(`CMDS ↑`)
  })

  it(`toggles discovery set when clicking buttons`, async () => {
    const [full_test_btn, unique_protos_btn] = toggle_buttons(`Full Test Set`)

    expect(unique_protos_btn.classList.contains(`active`)).toBe(true)
    expect(full_test_btn.classList.contains(`active`)).toBe(false)

    full_test_btn.click()
    await tick()

    expect(unique_protos_btn.classList.contains(`active`)).toBe(false)
    expect(full_test_btn.classList.contains(`active`)).toBe(true)
  })

  it(`updates column visibility when toggling checkboxes`, async () => {
    const columns_btn = doc_query<HTMLButtonElement>(`details.column-toggles summary`)
    const column_menu = doc_query(`.column-menu`)
    const details = doc_query<HTMLDetailsElement>(`details.column-toggles`)

    expect(column_menu.closest(`details`)).toBe(details)
    expect(details.open).toBe(false)

    columns_btn.click()
    expect(details.open).toBe(true)

    let f1_cells = document.querySelectorAll(`th, td`)
    let has_f1_column = [...f1_cells].some((cell) => cell.textContent?.includes(`F1`))
    expect(has_f1_column).toBe(true)

    const checkboxes = document.querySelectorAll<HTMLInputElement>(
      `.column-menu input[type="checkbox"]`,
    )
    const f1_checkbox = [...checkboxes].find((cb) =>
      cb.parentElement?.textContent?.includes(`F1`),
    )
    if (!f1_checkbox) throw new Error(`F1 checkbox not found`)
    expect(f1_checkbox.checked).toBe(true)

    f1_checkbox.click()
    await tick()
    expect(f1_checkbox.checked).toBe(false)

    f1_cells = document.querySelectorAll(`th, td`)
    has_f1_column = [...f1_cells].some((cell) => cell.textContent?.includes(`F1`))
    expect(has_f1_column).toBe(false)

    columns_btn.click()
    expect(details.open).toBe(false)
  })

  it(`toggles non-compliant models`, async () => {
    const selected_scatter_label = () =>
      doc_query(`span.selected-label`).textContent?.replaceAll(/\s+/g, ` `)
    const toggle = doc_query<HTMLInputElement>(
      `.table-controls label input[type="checkbox"]`,
    )
    expect(toggle.parentElement?.textContent).toContain(`Compliant models`)
    expect(toggle.checked).toBe(true)

    const n_models_on_load = document.querySelectorAll(`tbody tr`).length
    expect(selected_scatter_label()).toContain(`Params`)
    expect(selected_scatter_label()).toContain(`${n_models_on_load} models`)

    toggle.click()
    expect(toggle.checked).toBe(false)
    await tick()
    const n_non_compliant_only = document.querySelectorAll(`tbody tr`).length
    expect(n_non_compliant_only).toBeLessThan(n_models_on_load)
    expect(selected_scatter_label()).toContain(`${n_non_compliant_only} models`)
  })

  it(`displays valid best model information`, () => {
    const text = doc_query(`#best-report`).textContent ?? ``
    expect(text).toMatch(/highest F1 score/)
    expect(text).toMatch(/discovery acceleration factor/)

    const f1_match = /F1 score of (?<f1>[\d.]+)/.exec(text)
    const daf_match = /DAF\) of (?<daf>[\d.]+)/.exec(text)
    if (!f1_match?.groups || !daf_match?.groups) {
      throw new Error(`Could not find F1 or DAF values in text: ${text}`)
    }
    const f1_val = Number(f1_match.groups.f1)
    const daf_val = Number(daf_match.groups.daf)
    expect(f1_val > 0 && f1_val < 1, `F1=${f1_val} is out of range`).toBe(true)
    expect(daf_val > 0 && daf_val < 10, `DAF=${daf_val} is out of range`).toBe(true)
  })

  it(`renders table downloads section`, () => {
    // Check that the download buttons section exists
    const download_section = doc_query(`.downloads`)

    // Check that it contains SVG, PNG, CSV, Excel buttons plus RSS link
    const download_buttons = download_section.querySelectorAll(`.download-btn`)
    expect(download_buttons).toHaveLength(5)

    const buttons = [...download_buttons].map((btn) => btn.textContent?.trim())
    expect(buttons).toStrictEqual([`SVG`, `PNG`, `CSV`, `Excel`, `RSS`])
  })
})

describe(`Landing Page URL state`, () => {
  const column_toggle = (): HTMLInputElement => {
    const checkbox = document.querySelector<HTMLInputElement>(
      `details.column-toggles input:not(:disabled)`,
    )
    if (!checkbox) throw new Error(`No column toggle checkbox found`)
    return checkbox
  }
  const reset_columns_button = (): HTMLButtonElement => {
    const button = document.querySelector<HTMLButtonElement>(
      `button[aria-label="Reset all columns to defaults"]`,
    )
    if (!button) throw new Error(`No reset columns button found`)
    return button
  }

  it.each([
    [`http://localhost/?preset=MD&sort=F1`, `F1`, `descending`],
    [`http://localhost/?preset=MD&sort=combined_score&dir=asc`, `CMDS`, `ascending`],
  ])(
    `preserves URL sort when restoring a column preset`,
    async (url, sorted_col, aria_sort) => {
      await mount_with_url(Page, url)

      expect(header_text()).toContain(`vDOS`)
      expect(sorted_header()?.textContent).toContain(sorted_col)
      expect(sorted_header()?.getAttribute(`aria-sort`)).toBe(aria_sort)
    },
  )

  it(`omits preset default sort params from the URL`, async () => {
    await mount_with_url(Page, `http://localhost/?preset=MD&sort=combined_score&dir=desc`)

    expect(location.search).toBe(`?preset=MD`)
    expect(sorted_header()?.textContent).toContain(`CMDS`)
    expect(sorted_header()?.getAttribute(`aria-sort`)).toBe(`descending`)

    preset_button(`Geo Opt`).click()
    await tick()
    expect(sorted_header()?.textContent).toContain(`RMSD`)
    expect(sorted_header()?.getAttribute(`aria-sort`)).toBe(`ascending`)
  })

  it(`drops preset from the URL only while columns are customized`, async () => {
    await mount_with_url(Page, `http://localhost/?preset=MD`)
    expect(location.search).toContain(`preset=MD`)

    column_toggle().click()
    await tick()
    expect(location.search).not.toContain(`preset=`)

    reset_columns_button().click()
    await tick()
    expect(location.search).toContain(`preset=MD`)
  })
})
