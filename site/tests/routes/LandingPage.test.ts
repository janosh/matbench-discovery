import Page from '$routes/+page.svelte'
import { tick } from 'svelte'
import { beforeEach, describe, expect, it } from 'vitest'
import { doc_query, mount, mount_with_url, sorted_header } from '../index'

const header_text = () =>
  [...document.querySelectorAll(`thead th`)].map((header) => header.textContent).join(` `)
const toggle_buttons = (label: string): HTMLButtonElement[] => {
  const toggle = [...document.querySelectorAll(`.selection-toggle`)].find((element) =>
    element.textContent?.includes(label),
  )
  if (!toggle) throw new Error(`No selection toggle contains ${label}`)
  return [...toggle.querySelectorAll<HTMLButtonElement>(`button`)]
}
const preset_button = (label: string): HTMLButtonElement => {
  const button = toggle_buttons(label).find(
    (candidate) => candidate.textContent?.trim() === label,
  )
  if (!button) throw new Error(`No preset button found for ${label}`)
  return button
}
const table_header = (label: string): HTMLTableCellElement => {
  const header = [...document.querySelectorAll<HTMLTableCellElement>(`thead th`)].find(
    (candidate) => candidate.textContent?.replace(/\s*[↑↓]\s*$/, ``).trim() === label,
  )
  if (!header) throw new Error(`no table header labeled ${label}`)
  return header
}
const expect_sort = (label: string, direction: `ascending` | `descending`) => {
  const header = sorted_header()
  expect(header?.textContent).toContain(label)
  expect(header?.getAttribute(`aria-sort`)).toBe(direction)
}

describe(`Landing Page`, () => {
  // happy-dom mounts of the full-column metrics table are slow in CI
  beforeEach(() => {
    mount(Page, { target: document.body })
  }, 30_000)

  const select_preset = async (label: string) => {
    preset_button(label).click()
    await tick()
  }

  it(`renders Discovery metrics by default`, () => {
    const preset_labels = toggle_buttons(`Discovery`).map((button) =>
      button.textContent?.trim(),
    )
    expect(preset_labels).toStrictEqual([
      `Discovery`,
      `Geo Opt`,
      `Phonons`,
      `MD`,
      `Diatomics`,
    ])
    expect(document.body.textContent).toContain(`Discovery test set:`)
    expect(header_text()).toContain(`DAF`)
    expect(header_text()).toContain(`Params`)
    expect(header_text()).toContain(`rcut`)
    expect(doc_query(`#best-report`).textContent).toContain(`Discovery view`)
    expect(doc_query(`#best-report`).textContent).toContain(`F1`)
  })

  it.each([
    [`Geo Opt`, `Geometry Optimization`],
    [`MD`, `Molecular Dynamics`],
  ])(`expands %s in its tooltip`, async (label, expanded_label) => {
    const button = preset_button(label)
    button.dispatchEvent(new MouseEvent(`mouseenter`))
    await new Promise((resolve) => setTimeout(resolve, 110))

    const tooltip_id = button.getAttribute(`aria-describedby`)
    expect(tooltip_id).not.toBeNull()
    expect(document.querySelector(`[id="${tooltip_id}"]`)?.textContent).toContain(
      expanded_label,
    )
  })

  // Each non-default task preset reveals one of its signature columns.
  it.each([
    [`Phonons`, [`κSRE`, `κSRME`, `κSRD`, `κ failed`, `Im(ω)`, `W1(ω)`], `κSRME`], // all six phonon metrics
    [`Geo Opt`, [`Σ`], `RMSD`], // symmetry metrics (Σ= / Σ↓ / Σ↑)
    [`MD`, [`vDOS`], `CMDS`], // vDOS err (RDF is hidden from leaderboards as redundant)
    [`Diatomics`, [`E jump`], `CDS`],
  ])(
    `%s preset updates task metrics, headlines and best model`,
    async (preset, signature_columns, primary_metric) => {
      expect(header_text()).not.toContain(signature_columns[0])

      await select_preset(preset)
      const headers = header_text()
      const best_report = doc_query(`#best-report`).textContent
      for (const signature_column of signature_columns) {
        expect(headers).toContain(signature_column)
      }
      for (const headline of [`CPS`, `F1`, `RMSD`, `CMDS`, `CDS`]) {
        expect(headers).toContain(headline)
      }
      expect(best_report).toContain(`${preset} view`)
      expect(best_report).toContain(primary_metric)
    },
  )

  it(`surfaces a beta warning for the MD metrics`, async () => {
    const beta_warning_count = () =>
      [...document.querySelectorAll(`blockquote`)].filter((blockquote) =>
        blockquote.textContent?.toLowerCase().includes(`interpret with caution`),
      ).length
    expect(beta_warning_count()).toBe(1) // always shown in the page's MD note

    await select_preset(`MD`)
    expect(beta_warning_count()).toBe(2) // + contextual warning above the MD table
  })

  it(`toggles the discovery set and hides it outside Discovery`, async () => {
    const test_set_shown = () =>
      [...document.querySelectorAll(`.selection-toggle`)].some((toggle) =>
        toggle.textContent?.includes(`Full Test Set`),
      )
    const [full_test_button, unique_protos_button] = toggle_buttons(`Full Test Set`)
    expect(test_set_shown()).toBe(true)
    expect(unique_protos_button.classList.contains(`active`)).toBe(true)
    expect(full_test_button.classList.contains(`active`)).toBe(false)

    full_test_button.click()
    await tick()
    expect(unique_protos_button.classList.contains(`active`)).toBe(false)
    expect(full_test_button.classList.contains(`active`)).toBe(true)

    await select_preset(`MD`)
    expect(test_set_shown()).toBe(false)
  })

  it(`auto-sorts presets until the user manually sorts the table`, async () => {
    expect(header_text()).toContain(`F1 ↑`)

    doc_query(`section.full-bleed`).click()
    await tick()
    await select_preset(`MD`)
    expect(header_text()).toContain(`CMDS ↑`)

    table_header(`F1`).click()
    await tick()
    expect(header_text()).toMatch(/F1 [↑↓]/)

    await select_preset(`Geo Opt`)
    expect(header_text()).toMatch(/F1 [↑↓]/)
    expect(header_text()).not.toContain(`RMSD ↓`)
  })

  it(`updates column visibility when toggling checkboxes`, async () => {
    const columns_button = doc_query<HTMLButtonElement>(`details.column-toggles summary`)
    const column_menu = doc_query(`.column-menu`)
    const details = doc_query<HTMLDetailsElement>(`details.column-toggles`)

    expect(column_menu.parentElement).toBe(document.body)
    expect(details.open).toBe(false)

    columns_button.click()
    expect(details.open).toBe(true)

    expect(header_text()).toContain(`F1`)

    const checkboxes = document.querySelectorAll<HTMLInputElement>(
      `.column-menu input[type="checkbox"]`,
    )
    const f1_checkbox = [...checkboxes].find((checkbox) =>
      checkbox.parentElement?.textContent?.includes(`F1`),
    )
    if (!f1_checkbox) throw new Error(`F1 checkbox not found`)
    expect(f1_checkbox.checked).toBe(true)

    f1_checkbox.click()
    await tick()
    expect(f1_checkbox.checked).toBe(false)

    expect(header_text()).not.toContain(`F1`)

    columns_button.click()
    expect(details.open).toBe(false)
  })

  it(`filters models via the training-data dropdown`, async () => {
    const selected_scatter_label = () =>
      doc_query(`span.selected-label`).textContent?.replaceAll(/\s+/g, ` `)
    const model_count_on_load = document.querySelectorAll(`tbody tr`).length
    expect(selected_scatter_label()).toContain(`Params`)
    expect(selected_scatter_label()).toContain(`${model_count_on_load} models`)

    // Excluding OMat24 removes many otherwise visible models.
    const training_menu = [...document.querySelectorAll(`details.filter-menu`)].find(
      (menu) => menu.querySelector(`summary`)?.textContent?.includes(`Training`),
    )
    if (!training_menu) throw new Error(`Training data filter menu not found`)
    const omat_exclude_checkbox = training_menu.querySelector<HTMLInputElement>(
      `input[aria-label="exclude OMat24"]`,
    )
    if (!omat_exclude_checkbox) throw new Error(`OMat24 exclude checkbox not found`)
    omat_exclude_checkbox.click()
    await tick()

    const filtered_model_count = document.querySelectorAll(`tbody tr`).length
    expect(filtered_model_count).toBeLessThan(model_count_on_load)
    expect(selected_scatter_label()).toContain(`${filtered_model_count} models`)
  })

  it(`renders table downloads section`, () => {
    const download_section = doc_query(`.downloads`)
    const download_buttons = download_section.querySelectorAll(`.download-btn`)
    expect(download_buttons).toHaveLength(5)

    const buttons = [...download_buttons].map((button) => button.textContent?.trim())
    expect(buttons).toStrictEqual([`SVG`, `PNG`, `CSV`, `Excel`, `RSS`])
  })
})

describe(`Landing Page URL state`, () => {
  it.each([
    [`http://localhost/?preset=MD&sort=F1`, `F1`, `descending`],
    [`http://localhost/?preset=MD&sort=combined_score&dir=asc`, `CMDS`, `ascending`],
  ] as const)(
    `preserves URL sort when restoring a column preset`,
    async (url, sorted_column, aria_sort) => {
      await mount_with_url(Page, url)

      expect(header_text()).toContain(`vDOS`)
      expect_sort(sorted_column, aria_sort)
    },
  )

  it(`canonicalizes preset URL across default and customized columns`, async () => {
    await mount_with_url(Page, `http://localhost/?preset=MD&sort=combined_score&dir=desc`)

    expect(location.search).toBe(`?preset=MD`)
    expect_sort(`CMDS`, `descending`)

    preset_button(`Geo Opt`).click()
    await tick()
    expect_sort(`RMSD`, `ascending`)

    doc_query<HTMLInputElement>(
      `.column-menu input[type="checkbox"]:not(:disabled)`,
    ).click()
    await tick()
    expect(new URLSearchParams(location.search).has(`preset`)).toBe(false)

    doc_query<HTMLButtonElement>(
      `button[aria-label="Reset all columns to defaults"]`,
    ).click()
    await tick()
    expect(new URLSearchParams(location.search).get(`preset`)).toBe(`Geo Opt`)
  })

  it(`restores the heatmap toggle from the URL and omits it at default`, async () => {
    await mount_with_url(Page, `http://localhost/?heatmap=0`)

    const heatmap_toggle = doc_query<HTMLInputElement>(
      `input[aria-label="Toggle heatmap colors"]`,
    )
    expect(heatmap_toggle.checked).toBe(false)

    heatmap_toggle.click()
    await tick()
    expect(location.search).not.toContain(`heatmap=`)
  })
})
