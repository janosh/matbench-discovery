import { OPENNESS_OPTIONS } from '$lib/url-state.svelte'
import { comparison } from '$lib/model-comparison.svelte'
import Page from '$routes/+page.svelte'
import { flushSync, tick } from 'svelte'
import { beforeEach, describe, expect, it, vi } from 'vitest'
import {
  doc_query,
  header_name,
  mount,
  mount_with_url,
  navigate,
  POPOVER_OPEN_ATTR,
  query_param,
  sorted_header,
} from '../index'

const header_text = () =>
  [...document.querySelectorAll(`thead th`)]
    .map((header) => header.textContent)
    .join(` `)
    .replaceAll(/\s+/g, ` `)
const toggle_buttons = (label: string): HTMLButtonElement[] => {
  const toggle = [...document.querySelectorAll(`.button-group`)].find((element) =>
    element.textContent?.includes(label),
  )
  if (!toggle) throw new Error(`No button group contains ${label}`)
  return [...toggle.querySelectorAll<HTMLButtonElement>(`button`)]
}
// label of the active button in the toggle group containing `label`
const pressed_toggle = (label: string): string | undefined =>
  toggle_buttons(label)
    .find((button) => button.getAttribute(`aria-checked`) === `true`)
    ?.textContent?.trim()
const preset_button = (label: string): HTMLButtonElement => {
  const button = toggle_buttons(label).find(
    (candidate) => candidate.textContent?.trim() === label,
  )
  if (!button) throw new Error(`No preset button found for ${label}`)
  return button
}
const table_header = (label: string): HTMLTableCellElement => {
  const header = [...document.querySelectorAll<HTMLTableCellElement>(`thead th`)].find(
    (candidate) => header_name(candidate) === label,
  )
  if (!header) throw new Error(`no table header labeled ${label}`)
  return header
}
const expect_sort = (label: string, direction: `ascending` | `descending`) => {
  const header = sorted_header()
  expect(header?.textContent).toContain(label)
  expect(header?.getAttribute(`aria-sort`)).toBe(direction)
}

const mount_page = () => mount(Page, { target: document.body })

describe(`Landing Page`, () => {
  let page_component: ReturnType<typeof mount_page>
  // happy-dom mounts of the full-column metrics table are slow in CI
  beforeEach(() => {
    page_component = mount_page()
  }, 30_000)

  const select_preset = async (label: string) => {
    preset_button(label).click()
    await tick()
  }

  it(`keeps the full table, plots, score weights and About visible with table help collapsed`, async () => {
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
    expect(
      [...document.querySelectorAll(`.button-group [aria-label]`)].map((el) =>
        el.getAttribute(`aria-label`),
      ),
    ).toEqual(expect.arrayContaining([`Column presets`, `Discovery test set`]))
    expect([...document.querySelectorAll(`thead th`)].map(header_name)).toEqual([
      `#`,
      `Model`,
      `CPS`,
      `Acc`,
      `F1`,
      `DAF`,
      `Prec`,
      `MAE`,
      `R2`,
      `κSRME`,
      `RMSD`,
      `CMDS`,
      `CDS`,
      `Params`,
      `Targets`,
      `Date Added`,
      `Links`,
      `rcut`,
      `Training Set`,
      `Org`,
    ])
    expect(document.querySelectorAll(`div.scatter`)).toHaveLength(2)

    const plots = [...document.querySelectorAll(`.plot-section`)]
    expect(plots.map((section) => section.getAttribute(`aria-labelledby`))).toEqual([
      `cps-progress-over-time`,
      `github-activity`,
    ])
    expect(plots[0].nextElementSibling).toBe(plots[1])
    expect(doc_query(`p`, plots[1]).querySelector(`a`)).toBeNull()
    expect(document.querySelectorAll(`a[href="/benchmarks"]`)).toHaveLength(1)
    expect(document.querySelectorAll(`a[href="/contribute"]`)).toHaveLength(1)
    const details = [...document.querySelectorAll<HTMLDetailsElement>(`.page-details`)]
    expect(
      details.map((element) => element.querySelector(`summary`)?.textContent),
    ).toEqual([`How to read the table`])
    expect(details.map(({ open }) => open)).toEqual([false])
    expect(doc_query(`#score-weights`).closest(`details`)).toBeNull()
    expect(doc_query(`#score-weights-heading`).textContent).toBe(`Adjust score weights`)
    expect(doc_query(`#score-weights svg[aria-label^="Radar chart"]`)).toBeDefined()
    expect(doc_query(`#table-guide`).textContent).toContain(`missing-result explanations`)
    const about = doc_query(`#about-benchmark`)
    expect(about.closest(`details`)).toBeNull()
    expect(doc_query(`h2`, about).textContent).toBe(`About Matbench Discovery`)
    expect(
      [...about.querySelectorAll<HTMLAnchorElement>(`p:first-of-type a`)].map((link) => [
        link.textContent,
        link.getAttribute(`href`),
      ]),
    ).toEqual([
      [`crystal discovery`, `/benchmarks/discovery`],
      [`geometry optimization`, `/benchmarks/geo-opt`],
      [`phonons`, `/benchmarks/phonons`],
      [`molecular dynamics`, `/benchmarks/md`],
      [`diatomics`, `/benchmarks/diatomics`],
    ])
    expect(doc_query(`#about-benchmark a[href^="https://doi.org/"]`)).toBeDefined()
    expect(document.querySelector(`.table-footer`)).toBeNull()
    for (const [href, label] of [
      [`/contribute`, `Submit a model`],
      [`/rss.xml`, `RSS`],
    ]) {
      const link = doc_query(`.control-buttons a[href="${href}"]`)
      expect(link.textContent?.trim()).toBe(label)
      expect(link.compareDocumentPosition(doc_query(`table`))).toBe(
        Node.DOCUMENT_POSITION_FOLLOWING,
      )
    }
    expect(
      doc_query(`.control-buttons button[aria-haspopup="menu"]`).textContent?.trim(),
    ).toBe(`Export`)
  })

  it.each([
    [`Geo Opt`, `Geometry Optimization`],
    [`MD`, `Molecular Dynamics`],
  ])(`expands %s in its tooltip`, async (label, expanded_label) => {
    const button = preset_button(label)
    flushSync()
    button.dispatchEvent(new PointerEvent(`pointerover`))
    let tooltip: HTMLElement | undefined
    await vi.waitFor(() => {
      const tooltip_id = button.getAttribute(`aria-describedby`)
      expect(tooltip_id).not.toBeNull()
      tooltip = doc_query(`[id="${tooltip_id}"]`)
      expect(tooltip.textContent).toContain(expanded_label)
    })
    // svelte-widgets 1.6 shows the tooltip as a top-layer popover (polyfilled in tests/index.ts)
    expect(tooltip?.hasAttribute(`popover`)).toBe(true)
    expect(tooltip?.hasAttribute(POPOVER_OPEN_ATTR)).toBe(true)
    expect(tooltip?.hidden).toBe(false)

    button.dispatchEvent(new PointerEvent(`pointerout`))
    await vi.waitFor(() => {
      expect(button.hasAttribute(`aria-describedby`)).toBe(false)
      expect(tooltip?.hasAttribute(POPOVER_OPEN_ATTR)).toBe(false)
      expect(tooltip?.hidden).toBe(true)
    })
  })

  // Each non-default task preset reveals one of its signature columns.
  it.each([
    [`Phonons`, [`κSRE`, `κSRME`, `κSRD`, `κ failed`, `Im(ω)`, `W1(ω)`]], // all six phonon metrics
    [`Geo Opt`, [`Σ`]], // symmetry metrics (Σ= / Σ↓ / Σ↑)
    [`MD`, [`vDOS`]], // vDOS err (RDF is hidden from leaderboards as redundant)
    [`Diatomics`, [`E jump`]],
  ])(
    `%s preset reveals task metrics and keeps headline metrics`,
    async (preset, signature_columns) => {
      expect(header_text()).not.toContain(signature_columns[0])

      await select_preset(preset)
      const headers = header_text()
      for (const signature_column of signature_columns) {
        expect(headers).toContain(signature_column)
      }
      for (const headline of [`CPS`, `F1`, `RMSD`, `CMDS`, `CDS`]) {
        expect(headers).toContain(headline)
      }
    },
  )

  it(`shows a concise beta note only when the MD task is selected`, async () => {
    expect(document.querySelector(`.task-note`)).toBeNull()
    await select_preset(`MD`)
    const note = doc_query(`.task-note`)
    expect(note.textContent).toContain(`MD is in beta.`)
    expect(note.textContent).toContain(`preliminary and may change`)
    expect(doc_query(`a`, note).getAttribute(`href`)).toBe(`/benchmarks/md`)
    await select_preset(`Discovery`)
    expect(document.querySelector(`.task-note`)).toBeNull()
  })

  it(`toggles the discovery set and hides it outside Discovery`, async () => {
    expect(header_text()).toMatch(/F1.*DAF.*Prec/)
    expect(header_text()).not.toContain(`Recall`)

    const test_set_shown = () =>
      [...document.querySelectorAll(`.button-group`)].some((toggle) =>
        toggle.textContent?.includes(`Full Test Set`),
      )
    const [full_test_button] = toggle_buttons(`Full Test Set`)
    expect(test_set_shown()).toBe(true)
    expect(pressed_toggle(`Full Test Set`)).toBe(`Unique Prototypes`)

    full_test_button.click()
    await tick()
    expect(pressed_toggle(`Full Test Set`)).toBe(`Full Test Set`)

    await select_preset(`MD`)
    expect(test_set_shown()).toBe(false)
  })

  it(`ignores snapshot values that no longer exist`, async () => {
    const model_count = document.querySelectorAll(`tbody tr`).length

    // a snapshot captured before a deploy that renamed/removed these
    page_component.snapshot.restore({
      discovery_set: `most_stable_10k`,
      col_preset: `Removed Preset`,
      // openness left unrestricted so the row count only reflects the unknown dataset
      filters: {
        training: { RenamedDataset: `require` },
        openness: [...OPENNESS_OPTIONS],
      },
    })
    await tick()

    // the stale set and preset fall back to the defaults the page mounted with
    expect(pressed_toggle(`Full Test Set`)).toBe(`Unique Prototypes`)
    expect(pressed_toggle(`Discovery`)).toBe(`Discovery`)
    // the unknown dataset is ignored, leaving the row count unchanged
    expect(document.querySelectorAll(`tbody tr`)).toHaveLength(model_count)
  })

  it.each([`click`, `Enter`, ` `])(
    `auto-sorts presets until manual sorting with %j`,
    async (event_type) => {
      const cps_values = [...document.querySelectorAll(`td[data-col="CPS"]`)]
        .map((cell) => Number(cell.textContent?.trim()))
        .filter(Number.isFinite)
      expect(cps_values.length).toBeGreaterThan(1)
      expect(cps_values).toStrictEqual(
        [...cps_values].toSorted((value_1, value_2) => value_2 - value_1),
      )
      expect_sort(`CPS`, `descending`)

      doc_query(`section.full-bleed`).click()
      await tick()
      await select_preset(`MD`)
      expect(header_text()).toContain(`CMDS ↑`)

      const f1_header = table_header(`F1`)
      // Sortable headers retain their column-header semantics.
      f1_header.setAttribute(`role`, `columnheader`)
      if (event_type === `click`) f1_header.click()
      else
        f1_header.dispatchEvent(
          new KeyboardEvent(`keydown`, { key: event_type, bubbles: true }),
        )
      await tick()
      expect(header_text()).toMatch(/F1 [↑↓]/)

      await select_preset(`Geo Opt`)
      expect(header_text()).toMatch(/F1 [↑↓]/)
      expect(header_text()).not.toContain(`RMSD ↓`)

      preset_button(`Geo Opt`).focus()
      preset_button(`Geo Opt`).dispatchEvent(
        new KeyboardEvent(`keydown`, { key: `ArrowRight`, bubbles: true }),
      )
      await tick()
      expect(pressed_toggle(`Discovery`)).toBe(`Phonons`)
      expect(header_text()).toMatch(/F1 [↑↓]/)
    },
  )

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
    comparison.set([])
    const selected_scatter_label = () =>
      [...document.querySelectorAll(`.property-picker`)]
        .find((picker) => picker.querySelector(`label`)?.textContent === `Marker size`)
        ?.querySelector(`.selected-label`)
        ?.textContent?.replaceAll(/\s+/g, ` `)
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
    const model_count = doc_query(`.ranking-context strong`)
    expect(model_count.textContent).toMatch(
      new RegExp(`^${filtered_model_count} of \\d+ eligible models$`),
    )

    doc_query(`tbody tr`).dispatchEvent(new MouseEvent(`dblclick`, { bubbles: true }))
    await tick()
    const selected_only = doc_query<HTMLInputElement>(
      `input[aria-label="Toggle between showing only selected models and all models"]`,
    )
    selected_only.click()
    await tick()
    expect(document.querySelectorAll(`tbody tr`)).toHaveLength(1)
    expect(model_count.textContent).toMatch(/^1 of \d+ eligible models$/)
    comparison.set([])
    await tick()
    expect(model_count.textContent).toMatch(/^0 of \d+ eligible models$/)
    selected_only.click()
    await tick()
    expect(model_count.textContent).toMatch(
      new RegExp(`^${filtered_model_count} of \\d+ eligible models$`),
    )
  })
})

describe(`Landing Page URL state`, () => {
  const column_checkbox = () =>
    doc_query<HTMLInputElement>(`.column-menu input[type="checkbox"]:not(:disabled)`)

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

  it.each([
    [`Discovery`, `CPS`, null, `descending`],
    [`Geo Opt`, `RMSD`, `rmsd`, `ascending`],
    [`Phonons`, `κSRME`, `κ_SRME`, `ascending`],
    [`MD`, `CMDS`, `combined_score`, `descending`],
    [`Diatomics`, `CDS`, `diatomics_combined_score`, `descending`],
  ] as const)(
    `canonicalizes %s URLs and preserves sorting after customizing columns`,
    async (preset, sorted_column, _sort_key, direction) => {
      await mount_with_url(
        Page,
        `http://localhost/?preset=MD&sort=combined_score&dir=desc`,
      )

      expect(location.search).toBe(`?preset=MD`)
      expect_sort(`CMDS`, `descending`)

      preset_button(preset).click()
      await tick()
      expect_sort(sorted_column, direction)

      column_checkbox().click()
      await tick()
      const custom_query = location.search
      const custom_headers = [...document.querySelectorAll(`thead th`)].map(header_name)
      const cohort_size = document.querySelectorAll(`tbody tr`).length
      expect(query_param(`preset`)).toBe(preset === `Discovery` ? null : preset)
      expect(query_param(`columns`)).not.toBeNull()
      expect(query_param(`sort`)).toBeNull()
      expect(query_param(`dir`)).toBeNull()

      doc_query<HTMLButtonElement>(
        `button[aria-label="Reset all columns to defaults"]`,
      ).click()
      await tick()
      expect(query_param(`preset`)).toBe(preset === `Discovery` ? null : preset)
      expect(query_param(`sort`)).toBeNull()
      expect(query_param(`dir`)).toBeNull()

      await navigate(custom_query)
      expect(pressed_toggle(`Discovery`)).toBe(preset)
      expect_sort(sorted_column, direction)
      expect([...document.querySelectorAll(`thead th`)].map(header_name)).toEqual(
        custom_headers,
      )
      expect(document.querySelectorAll(`tbody tr`)).toHaveLength(cohort_size)
    },
  )

  it(`renders unknown URL sort labels as text`, async () => {
    const unknown_metric = `<b>unknown metric</b>`
    await mount_with_url(
      Page,
      `http://localhost/?sort=${encodeURIComponent(unknown_metric)}&columns=Model,F1`,
    )
    const context = doc_query(`.ranking-context`)
    expect(context.textContent).toContain(unknown_metric)
    expect(context.querySelector(`b`)).toBeNull()
  })

  it(`restores columns, filters, heatmap, and scores on initial and later navigation`, async () => {
    await mount_with_url(Page, `http://localhost/?heatmap=0`)

    const first_row = doc_query(`tbody tr`)
    const cps_cell = doc_query(`td[data-col="CPS"]`, first_row)
    const heatmap_toggle = doc_query<HTMLInputElement>(
      `input[aria-label="Toggle heatmap colors"]`,
    )
    expect(heatmap_toggle.checked).toBe(false)
    expect(cps_cell.style.getPropertyValue(`--cell-bg`)).toBe(``)

    heatmap_toggle.click()
    await tick()
    expect(location.search).not.toContain(`heatmap=`)
    expect(cps_cell.style.getPropertyValue(`--cell-bg`)).not.toBe(``)

    const filtered_query = `?train=MPtrj&openness=OSOD&targets=S&set=full_test_set&preset=MD&heatmap=0`
    // Revisit the same preset as well as a different one after customizing columns.
    for (const query of [filtered_query, filtered_query, ``]) {
      column_checkbox().click()
      await tick()
      await navigate(query, `popstate`)
      expect(location.search).toBe(query)
      expect(column_checkbox().checked).toBe(true)
      expect(header_text()).toContain(`Model`)
      if (query) expect(header_text()).toContain(`vDOS`)
      expect(heatmap_toggle.checked).toBe(!query)
      expect(
        doc_query<HTMLInputElement>(`.desktop-filters [aria-label="require MPtrj"]`)
          .checked,
      ).toBe(Boolean(query))
      expect(
        doc_query<HTMLInputElement>(`.desktop-filters [aria-label="require stress"]`)
          .checked,
      ).toBe(Boolean(query))
      expect(pressed_toggle(`Discovery`)).toBe(query ? `MD` : `Discovery`)
      expect_sort(query ? `CMDS` : `CPS`, `descending`)
    }
  })
})
