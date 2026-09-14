import { goto } from '$app/navigation'
import { page } from '$app/state'
import data_files from '$pkg/data-files.yml'
import { DIATOMICS_METRICS, DISCOVERY_SET_LABELS, HYPERPARAMS } from '$lib/labels'
import { comparison } from '$lib/model-comparison.svelte'
import {
  ACTIVE_MODELS,
  get_pred_file_urls,
  make_table_filters,
  score_weight_records,
} from '$lib/models.svelte'
import MetricsTable from '$lib/table/MetricsTable.svelte'
import type { DiscoverySet, Label, ModelData } from '$lib/types'
import { tick, type ComponentProps } from 'svelte'
import { beforeEach, describe, expect, it, onTestFinished, vi } from 'vitest'
import {
  doc_query,
  header_name,
  mount,
  mount_with_url,
  navigate,
  query_param,
} from '../index'
import app_css from '$site/src/app.css?raw'

const mount_table = (props: ComponentProps<typeof MetricsTable> = {}) =>
  mount(MetricsTable, { target: document.body, props })

// all header cells except the structural rank (#) column HeatmapTable renders
// for show_row_numbers
const header_cells = () => [
  ...document.querySelectorAll<HTMLTableCellElement>(`th:not(.row-num-col)`),
]
const header_names = () => header_cells().map(header_name)
const ranking_context = () =>
  doc_query(`.ranking-context`).textContent?.replaceAll(/\s+/g, ` `)
const column_checkbox = (name: string): HTMLInputElement => {
  const label = [...document.querySelectorAll(`.column-menu label`)].find(
    (candidate) => candidate.textContent?.trim() === name,
  )
  if (!label) throw new Error(`No column checkbox labeled ${name}`)
  return doc_query<HTMLInputElement>(`input`, label)
}

// table filters restricted to models trained (at least in part) on MPtrj
const mptrj_only_filters = () => {
  const filters = make_table_filters()
  filters.training = { MPtrj: `require` }
  return filters
}

const visible_row_count = (
  extra_filter: (model: ModelData) => boolean = make_table_filters().matches,
) => ACTIVE_MODELS.filter(extra_filter).length

// table filters with the default require-forces constraint cleared (shows all models
// incl. energy-only ones)
const all_targets_filters = () => {
  const filters = make_table_filters()
  filters.targets = {}
  return filters
}

describe(`MetricsTable`, () => {
  // the module-level comparison store persists across tests
  beforeEach(() => {
    comparison.keys.clear()
    comparison.open = false
  })
  // model behind each rendered row, resolved from its Model-cell link
  const row_models = (): ModelData[] =>
    [...document.querySelectorAll(`tbody td[data-col="Model"]`)].map((cell) => {
      const link = doc_query(`a[href^="/models/"]`, cell)
      const model_key = link.getAttribute(`href`)?.slice(`/models/`.length)
      const model = ACTIVE_MODELS.find((md) => md.model_key === model_key)
      if (!model) throw new Error(`no model for row link ${link.outerHTML}`)
      return model
    })

  it(
    `renders columns, header tooltips and prediction downloads`,
    { timeout: 30_000 },
    async () => {
      let discovery_set = $state<DiscoverySet>(`unique_prototypes`)
      mount_table({
        col_filter: () => true,
        get discovery_set() {
          return discovery_set
        },
      })

      const table = doc_query(`table`)
      expect(table.querySelector(`thead`)).not.toBeNull()
      expect(table.querySelector(`tbody`)).not.toBeNull()
      const table_container = doc_query(`.table-container`)
      expect(table_container.style.getPropertyValue(`--heatmap-column-max-width`)).toBe(
        `14.4em`,
      )
      expect(
        table_container.style.getPropertyValue(`--heatmap-sticky-cell-odd-bg`),
      ).toContain(`linear-gradient`)
      expect(
        table_container.style.getPropertyValue(`--heatmap-row-num-padding-left`),
      ).toBe(`0`)

      const header_texts = header_cells().map((header) =>
        header.textContent?.trim().replaceAll(/\s+/g, ` `),
      )
      // Model stays first and Org is a regular metadata column at the far right.
      expect(header_texts[0]).toBe(`Model`)
      expect(header_texts.at(-1)).toBe(`Org`)
      const metric_order = [`CPS ↑`, `F1`, `DAF`].map((col) => header_texts.indexOf(col))
      expect(metric_order).toStrictEqual([...metric_order].toSorted((n1, n2) => n1 - n2))

      const initial_rows = [...table.querySelectorAll(`tbody tr`)]
      expect(ranking_context()).toContain(`Sorted by CPS (descending)`)
      expect(ranking_context()).toContain(`MD and diatomics are excluded`)
      expect(ranking_context()).toContain(`CPS uses unique-prototype discovery scores`)
      for (const subset of [`full_test_set`, `unique_prototypes`] as const) {
        discovery_set = subset
        await tick()
        expect(ranking_context()).toContain(
          `Discovery: ${DISCOVERY_SET_LABELS[subset].label}`,
        )
        const models = row_models()
        for (const [idx, row] of [...table.querySelectorAll(`tbody tr`)].entries()) {
          expect(row).toBe(initial_rows[idx])
          const f1 = models[idx].metrics?.discovery?.[subset]?.F1
          expect(
            doc_query(`td[data-col="F1"]`, row).getAttribute(`data-sort-value`),
          ).toBe(f1 == null ? null : String(f1))
        }
      }

      const org_cell = doc_query(`td[data-col="Org"]:has(.org-preview)`)
      expect(org_cell.getAttribute(`style`)).not.toContain(`min-width:`)
      const cps_header = header_cells().find((header) => header_name(header) === `CPS`)
      if (!cps_header) throw new Error(`CPS header is missing`)
      const style = document.createElement(`style`)
      // happy-dom's CSS parser swallows the first rule after a bare @import.
      style.textContent = app_css.replaceAll(/^@import[^;]+;/gm, ``)
      document.body.append(style)
      cps_header.style.fontWeight = `700`
      await tick()
      expect(getComputedStyle(doc_query(`.control-buttons`)).alignItems).toBe(`baseline`)
      expect(doc_query(`.control-buttons`).firstElementChild).toBe(
        doc_query(`.control-buttons > [aria-label="Active model filters"]`),
      )
      for (const action of document.querySelectorAll(
        `.control-buttons > :is(button, a)`,
      )) {
        expect(getComputedStyle(action).whiteSpace).toBe(`nowrap`)
      }
      const trigger = doc_query<HTMLButtonElement>(`button[aria-haspopup]`, cps_header)
      trigger.dispatchEvent(new MouseEvent(`mouseenter`))
      await vi.waitFor(() => {
        const content = doc_query(`.popover`, cps_header)
        expect(content.textContent).toContain(`Combined Performance Score`)
        expect(content.textContent).toContain(`(higher=better)`)
        expect(getComputedStyle(content).fontWeight).toBe(`400`)
        expect(getComputedStyle(content).whiteSpace).toBe(`normal`)
        expect(getComputedStyle(content).textAlign).toBe(`left`)
      })
      trigger.dispatchEvent(new MouseEvent(`mouseleave`))
      await vi.waitFor(() => expect(cps_header.querySelector(`.popover`)).toBeNull())

      const pred_files_button = doc_query<HTMLButtonElement>(
        `tbody button[aria-label="Download model prediction files"]`,
      )
      expect(document.querySelector(`.pred-files-dropdown`)).toBeNull()

      // Exercise reopening and both dismissal paths; click_outside uses pointerdown.
      for (const [target, event] of [
        [document.body, new PointerEvent(`pointerdown`, { bubbles: true })],
        [globalThis, new KeyboardEvent(`keydown`, { key: `Escape` })],
      ] as const) {
        pred_files_button.click()
        await tick()
        expect(doc_query(`.pred-files-dropdown`).textContent).toContain(`Files for`)
        target.dispatchEvent(event)
        await tick()
        expect(document.querySelector(`.pred-files-dropdown`)).toBeNull()
      }
    },
  )

  it(`explains every n/a cell and displays the reason on hover`, async () => {
    const pending_model = ACTIVE_MODELS.find(
      ({ model_key }) => model_key === `equflash-29m-oam`,
    )
    if (!pending_model?.metrics) throw new Error(`Missing EquFlash test metrics`)
    const { metrics } = pending_model
    const previous_md = metrics.md
    metrics.md = undefined
    onTestFinished(() => {
      metrics.md = previous_md
    })
    mount_table({
      model_filter: (model: ModelData) =>
        [`cgcnn`, `equflash-29m-oam`].includes(model.model_key),
      filters: all_targets_filters(),
      col_filter: () => true,
    })
    await tick()
    const missing_cells = [...document.querySelectorAll(`tbody td`)].filter(
      (cell) => cell.textContent?.trim() === `n/a`,
    )
    expect(missing_cells.length).toBeGreaterThan(0)
    for (const cell of missing_cells) {
      const trigger = doc_query(`[data-title]`, cell)
      expect(trigger.getAttribute(`data-title`)).not.toBe(`Not available`)
      expect(trigger.getAttribute(`data-title`)?.length).toBeGreaterThan(15)
    }
    const unsupported = doc_query(`a[href="/models/cgcnn"]`).closest(`tr`)
    const pending = doc_query(`a[href="/models/equflash-29m-oam"]`).closest(`tr`)
    if (!unsupported || !pending) throw new Error(`Missing model rows`)
    expect(
      doc_query(`td[data-col="RMSD"] [data-title]`, unsupported).dataset.title,
    ).toContain(`requires forces`)
    const pending_cell = doc_query(`td[data-col="CMDS"]`, pending)
    expect(pending_cell.getAttribute(`data-sort-value`)).toBeNull()
    const trigger = doc_query(`[data-title]`, pending_cell)
    trigger.dispatchEvent(new PointerEvent(`pointerover`, { bubbles: true }))
    await vi.waitFor(() => {
      expect(doc_query(`.custom-tooltip .tooltip-content`).textContent).toBe(
        `Molecular Dynamics: not evaluated yet. Contributions welcome to add missing model predictions.`,
      )
    })
  })

  it.each([
    {
      n_valid: 66,
      n_eligible: 73,
      failed_elements: [`Cs`, `K`, `Rb`],
      missing_elements: [`Ac`, `Pa`, `Th`, `U`],
      tooltip: `Valid fits: 66/73 reference-eligible elements. No valid fit: 3 (Cs, K, Rb). Unavailable curves: 4 (Ac, Pa, Th, U).`,
    },
    {
      n_valid: 0,
      n_eligible: 2,
      failed_elements: [`K`],
      missing_elements: [`U`],
      tooltip: `Valid fits: 0/2 reference-eligible elements. No valid fit: 1 (K). Unavailable curves: 1 (U).`,
    },
    {
      n_valid: 2,
      n_eligible: 2,
      failed_elements: [],
      missing_elements: [],
      tooltip: `Valid fits: 2/2 reference-eligible elements. No valid fit: 0. Unavailable curves: 0.`,
    },
  ])(
    `shows frequency fit coverage $n_valid/$n_eligible without changing numeric sorting`,
    async ({ tooltip, ...coverage }) => {
      const model = ACTIVE_MODELS[0]
      const { metrics } = model
      if (!metrics) throw new Error(`Missing test model metrics`)
      const previous = metrics.diatomics
      onTestFinished(() => {
        metrics.diatomics = previous
      })
      metrics.diatomics = {
        pbe_vib_freq_error: coverage.n_valid ? 47.7 : undefined,
        pbe_vib_freq_coverage: coverage,
      }
      mount_table({
        model_filter: (candidate: ModelData) => candidate.model_key === model.model_key,
        column_labels: [DIATOMICS_METRICS.pbe_vib_freq_error],
        filters: all_targets_filters(),
      })
      await tick()
      const cell = doc_query(`tbody td[data-col="PBE Δω"]`)
      expect(cell.getAttribute(`data-sort-value`)).toBe(coverage.n_valid ? `47.7` : null)
      expect(cell.textContent?.trim()).toBe(
        `${coverage.n_valid ? `47.7` : `n/a`} · ${coverage.n_valid}/${coverage.n_eligible}`,
      )
      expect(doc_query(`small`, cell).style.fontWeight).toBe(`400`)
      doc_query(`[data-title]`, cell).dispatchEvent(
        new PointerEvent(`pointerover`, { bubbles: true }),
      )
      await vi.waitFor(() => {
        expect(doc_query(`.custom-tooltip .tooltip-content`).textContent).toBe(tooltip)
      })
    },
  )

  it.each([
    {
      name: `metadata`,
      hidden_keys: [`Training Set`, `Targets`, `benchmark_added`, `Links`],
      hidden_labels: [`Training Set`, `Targets`, `Date Added`, `Links`],
      retained_labels: [`CPS`, `F1`, `DAF`, `Prec`, `Acc`],
    },
    {
      name: `metrics`,
      hidden_keys: [`F1`, `DAF`],
      hidden_labels: [`F1`, `DAF`],
      retained_labels: [`Model`, `CPS`, `Prec`, `Acc`],
    },
  ])(
    `hides selected $name columns`,
    ({ hidden_keys, hidden_labels, retained_labels }) => {
      mount_table({
        col_filter: (col: Label) => !hidden_keys.includes(col.key ?? col.label),
      })

      const header_texts = header_names()
      for (const label of hidden_labels) expect(header_texts).not.toContain(label)
      for (const label of retained_labels) expect(header_texts).toContain(label)
    },
  )

  it(
    `hides energy-only models by default via the targets filter`,
    {
      timeout: 30_000,
    },
    async () => {
      // default filters require force prediction, hiding energy-only models
      const filters = make_table_filters()
      mount_table({ filters })
      await tick()
      const rows_without_energy = row_models().length
      expect(rows_without_energy).toBe(visible_row_count())

      // clearing the targets filter shows them
      filters.targets = {}
      await tick()
      const rows_with_energy = row_models().length

      expect(rows_with_energy).toBe(visible_row_count(filters.matches))
      expect(rows_with_energy).toBeGreaterThan(rows_without_energy)
    },
  )

  it(`reactively filters models and displays an empty state when none match`, async () => {
    let model_filter = $state<(model: ModelData) => boolean>(() => false)
    mount_table({
      get model_filter() {
        return model_filter
      },
    })
    const default_matches = make_table_filters().matches
    for (const matches of [
      () => false,
      () => true,
      (model: ModelData) => model.model_name.includes(`CHG`),
      () => false,
    ]) {
      model_filter = matches
      await tick()
      const expected_models = ACTIVE_MODELS.filter(
        (model) => default_matches(model) && matches(model),
      )
      expect(ranking_context()).toContain(
        `${expected_models.length} of ${ACTIVE_MODELS.filter(matches).length} eligible models`,
      )
      expect(
        row_models()
          .map(({ model_key }) => model_key)
          .toSorted(),
      ).toStrictEqual(expected_models.map(({ model_key }) => model_key).toSorted())
      expect(document.querySelectorAll(`tbody tr.empty-row`)).toHaveLength(
        expected_models.length ? 0 : 1,
      )
      if (!expected_models.length)
        expect(doc_query(`tbody`).textContent?.trim()).toBe(`No data`)
    }
  })

  it.each([
    {
      name: `sticky columns only`,
      col_filter: (col: Label) => col.sticky === true,
      expected_headers: [`Model`],
    },
    {
      name: `specific columns`,
      col_filter: (col: Label) => [`Model`, `F1`].includes(col.key ?? col.label),
      expected_headers: [`Model`, `F1`],
    },
    {
      name: `Model always first`,
      col_filter: (col: Label) => [`F1`, `Model`, `DAF`].includes(col.key ?? col.label),
      expected_headers: [`Model`, `F1`, `DAF`],
    },
  ])(`handles col_filter: $name`, ({ col_filter, expected_headers }) => {
    mount_table({ col_filter })

    expect(header_names()).toStrictEqual(expected_headers)
  })

  it.each([
    {
      model_filter: (model: ModelData) => model.model_name.includes(`CHG`),
      col_filter: (col: Label) => [`Model`, `F1`].includes(col.key ?? col.label),
      expected_model_match: `CHG`,
      expected_headers: [`Model`, `F1`],
    },
    {
      model_filter: (model: ModelData) => model.model_name.includes(`MACE`),
      col_filter: (col: Label) => [`Model`, `DAF`].includes(col.key ?? col.label),
      expected_model_match: `MACE`,
      expected_headers: [`Model`, `DAF`],
    },
  ])(
    `combines filters: $expected_model_match models with $expected_headers`,
    ({ model_filter, col_filter, expected_model_match, expected_headers }) => {
      mount_table({ model_filter, col_filter })

      expect(header_names()).toStrictEqual(expected_headers)

      const rows = document.querySelectorAll(`tbody tr`)
      const default_matches = make_table_filters().matches
      expect(rows).toHaveLength(
        visible_row_count((model) => default_matches(model) && model_filter(model)),
      )
      rows.forEach((row) => {
        const model_cell = row.querySelector(`td[data-col="Model"]`)
        expect(model_cell?.textContent).toContain(expected_model_match)
      })
    },
  )

  it(`marks models with excluded metric samples`, async () => {
    mount_table({
      model_filter: (model: ModelData) => model.model_name === `AlphaNet-v1-OAM`,
      col_filter: (col: Label) => col.label === `Model`,
    })
    await tick()

    const model_cell = doc_query(`td[data-col="Model"]`)
    const marker = model_cell.querySelector(
      `span[title="Diatomics metrics exclude He-He due to exploding errors"]`,
    )
    expect(model_cell.textContent).toContain(`AlphaNet-v1-OAM`)
    expect(marker?.textContent).toBe(`*`)
  })

  describe(`Column Sorting`, () => {
    // Date Added sorts by timestamp (chronological, not alphabetical); Training Set
    // and Params sort by their numeric data-sort-value, not display text
    it.each([
      {
        col_key: `benchmark_added`,
        header: `Date Added`,
        sort_key: (model: ModelData) => Date.parse(model.dates.benchmark_added ?? ``),
      },
      {
        col_key: `Training Set`,
        header: `Training Set`,
        sort_key: ({ n_training_materials, n_training_structures }: ModelData) =>
          n_training_materials == null || n_training_structures == null
            ? null
            : n_training_materials,
      },
      {
        col_key: HYPERPARAMS.model_params.key,
        header: `Params`,
        sort_key: (model: ModelData) => model.model_params,
      },
    ])(
      `sorts $header numerically via data-sort-value`,
      async ({ col_key, header, sort_key }) => {
        mount_table({
          col_filter: (col: Label) => [`Model`, col_key].includes(col.key ?? col.label),
        })
        await tick()

        const sort_header = header_cells().find((th) => th.textContent?.includes(header))
        if (!sort_header) throw new Error(`${header} column not found`)

        const cell_values = () =>
          row_models()
            .map(sort_key)
            .filter((value) => value != null)

        sort_header.click()
        await tick()
        expect(ranking_context()).toContain(`Sorted by ${header}`)
        expect(ranking_context()).not.toContain(`Discovery:`)
        expect(document.querySelector(`.cps-context`)).toBeNull()

        const values = cell_values()
        expect(values.length).toBeGreaterThan(1)
        // Known counts must be numerically sorted in either direction.
        const ascending = [...values].toSorted((val_1, val_2) => val_1 - val_2)
        expect([ascending, ascending.toReversed()]).toContainEqual(values)

        // second click toggles sort direction
        sort_header.click()
        await tick()
        expect(cell_values()).toStrictEqual(values.toReversed())
      },
    )

    it(`renders Training Set cells as dataset links`, async () => {
      mount_table({
        col_filter: (col: Label) => [`Model`, `Training Set`].includes(col.label),
      })
      await tick()

      const training_set_cells = [
        ...document.querySelectorAll(`td[data-col="Training Set"]`),
      ]
      expect(training_set_cells.length).toBeGreaterThan(0)
      for (const cell of training_set_cells) {
        // HTML cells get no td[data-sort-value] (HeatmapTable reads it off the inner
        // span instead), so a value there would mean raw markup leaked onto the td
        expect(cell.hasAttribute(`data-sort-value`)).toBe(false)
        expect(cell.querySelector(`a[href^="/data/"]`)).not.toBeNull()
      }
      const clipped = doc_query(
        `tr:has(a[href="/models/prophet-oame-mbd"]) td[data-col="Training Set"] .middle-ellipsis-html`,
      )
      expect(clipped.textContent).toBe(`16.4M (324M) MPtrj+OMat24+sAlex+ELEMENTA`)
      expect(clipped.lastElementChild?.textContent).toBe(`ELEMENTA`)
      expect(clipped.querySelector(`a[href="/data/elementa"]`)).not.toBeNull()
    })

    it.each([
      {
        test_name: `with all models shown`,
        props: { filters: all_targets_filters() },
      },
      {
        test_name: `with filtered columns`,
        props: {
          col_filter: (col: Label) => [`Model`, `CPS`, `F1`].includes(col.label),
        },
      },
      {
        test_name: `with an MPtrj-only training filter`,
        props: { filters: mptrj_only_filters() },
      },
    ])(
      `alphabetically sorts by Model name on $test_name header click`,
      { timeout: 30_000 }, // happy-dom renders of the full-column table are slow in CI
      async ({ props }) => {
        mount_table(props)
        await tick()

        const headers = header_cells()
        const model_header = headers.find((h) => h.textContent?.includes(`Model`))

        if (!model_header) throw new Error(`Model column header not found`)

        const get_model_names = () =>
          [...document.querySelectorAll(`td[data-col="Model"]`)]
            .map((cell) => {
              const link = cell.querySelector(`a`)
              return link?.getAttribute(`data-sort-value`)
            })
            .filter(Boolean) as string[]

        model_header.click() // sort ascending A-Z
        await tick()

        const sorted_model_names = get_model_names()

        expect(sorted_model_names).toHaveLength(
          visible_row_count(
            `filters` in props && props.filters ? props.filters.matches : undefined,
          ),
        )

        // alphabetical in either direction
        const ascending = sorted_model_names.toSorted((name_a, name_b) =>
          name_a.localeCompare(name_b),
        )
        expect([ascending, ascending.toReversed()]).toContainEqual(sorted_model_names)

        model_header.click()
        await tick()

        const reverse_sorted_model_names = get_model_names()
        // Second click should reverse the sort direction
        expect(reverse_sorted_model_names).toStrictEqual(sorted_model_names.toReversed())
      },
    )

    it(`prevents sorting of unsortable Links column`, async () => {
      mount_table({
        col_filter: (col: Label) =>
          [`Model`, `CPS`, `Links`].includes(col.key ?? col.label),
      })

      const headers = header_cells()
      const cps_header = headers.find((h) => h.textContent?.includes(`CPS`))
      if (!cps_header) throw new Error(`CPS column not found`)
      const links_header = headers.find((h) => h.textContent?.includes(`Links`))
      if (!links_header) throw new Error(`Links column not found`)

      expect(links_header.classList.contains(`not-sortable`)).toBe(true)

      // Sort by CPS first (to establish a known order)
      cps_header.click()
      await tick()

      const initial_models = [...document.querySelectorAll(`td[data-col="Model"]`)].map(
        (cell) => cell.textContent,
      )

      // Try to sort by Links
      links_header.click()
      await tick()

      const after_links_click_models = [
        ...document.querySelectorAll(`td[data-col="Model"]`),
      ].map((cell) => cell.textContent)

      // Order should not change
      expect(after_links_click_models).toStrictEqual(initial_models)
    })
  })

  describe(`Links Column`, () => {
    it(`renders external links, unavailable icons, and prediction file buttons`, async () => {
      const col_filter = (col: Label) => [`Model`, `Links`].includes(col.label)
      mount_table({ col_filter })

      await tick() // Wait for component to process data

      const links_cells = [...document.querySelectorAll(`td[data-col="Links"]`)]
      expect(links_cells).toHaveLength(visible_row_count())

      let rows_with_links = 0
      for (const cell of links_cells) {
        const links = [...cell.querySelectorAll(`a`)]
        if (links.length > 1) rows_with_links++

        for (const link of links) {
          expect(link.getAttribute(`target`)).toBe(`_blank`)
          expect(link.getAttribute(`rel`)).toBe(`noopener noreferrer`)

          const title = link.getAttribute(`title`)
          const href = link.getAttribute(`href`)
          expect(title).not.toBeNull()
          expect(title).not.toBe(``)
          expect(href).not.toBeNull()
          expect(href).not.toBe(``)

          // Each link should have a resolved SVG glyph, not an obsolete string icon.
          const svg = link.querySelector(`svg`)
          expect(svg?.getAttribute(`viewBox`)).toMatch(/\S/)
        }
      }

      // At least half of rows should have multiple links
      expect(rows_with_links).toBeGreaterThan(links_cells.length / 2)

      const missing_icon_titles = links_cells.flatMap((cell) =>
        [...cell.querySelectorAll(`span[title$="not available"] svg`)].map((icon) =>
          icon.closest(`span`)?.getAttribute(`title`),
        ),
      )

      // Note: this might fail if all models have all links, which is unlikely
      expect(missing_icon_titles.length).toBeGreaterThan(0)
      expect(missing_icon_titles.every((title) => title?.match(/not available/))).toBe(
        true,
      )

      // Find all pred_files buttons (every row renders one)
      const pred_file_buttons = [
        ...document.querySelectorAll(
          `button[aria-label="Download model prediction files"]`,
        ),
      ]
      expect(pred_file_buttons).toHaveLength(visible_row_count())

      for (const button of pred_file_buttons) {
        expect(button.getAttribute(`aria-label`)).toBe(`Download model prediction files`)

        const svg = button.querySelector(`svg`)
        expect(svg).not.toBeNull()
      }
    })
  })

  it(`renders the correct default columns`, () => {
    mount_table()

    // Core text expected in default visible columns (duplicates intended: MD and
    // diatomics each have Speed and Slowdown columns, disambiguated by tooltip)
    const expected_core_columns = [
      `Model`, // METADATA_COLS
      `Training Set`, // METADATA_COLS
      `Targets`, // METADATA_COLS
      `Date Added`, // METADATA_COLS
      `Links`, // METADATA_COLS
      `Org`, // METADATA_COLS
      `Params`, // HYPERPARAMS (short label)
      `rcut`, // HYPERPARAMS (short label) - textContent doesn't keep subscript
      `Σ= 10-2`, // ALL_METRICS (Geo Opt) - textContent doesn't keep superscript
      `Σ= 10-5`, // ALL_METRICS (Geo Opt) - textContent doesn't keep superscript
      `Σ↑ 10-2`, // ALL_METRICS (Geo Opt) - textContent doesn't keep superscript
      `Σ↑ 10-5`, // ALL_METRICS (Geo Opt) - textContent doesn't keep superscript
      `Σ↓ 10-2`, // ALL_METRICS (Geo Opt) - textContent doesn't keep superscript
      `Σ↓ 10-5`, // ALL_METRICS (Geo Opt) - textContent doesn't keep superscript
      `Acc`, // ALL_METRICS (Discovery - short)
      `F1`, // ALL_METRICS (Discovery - short)
      `DAF`, // ALL_METRICS (Discovery)
      `Prec`, // ALL_METRICS (Discovery - short)
      `Recall`, // ALL_METRICS (Discovery)
      `MAE`, // ALL_METRICS (Discovery)
      `R2`, // ALL_METRICS (Discovery)
      `RMSE`, // ALL_METRICS (Discovery)
      `κSRME`, // ALL_METRICS (Phonon) - textContent doesn't keep subscript
      `κSRE`, // ALL_METRICS (Phonon) - textContent doesn't keep subscript
      `κSRD`, // ALL_METRICS (Phonon) - textContent doesn't keep subscript
      `κ failed`, // ALL_METRICS (Phonon)
      `Im(ω)`, // ALL_METRICS (Phonon)
      `W1(ω)`, // ALL_METRICS (Phonon)
      `RMSD`, // ALL_METRICS (Geo Opt)
      `ΔERMSE`, // ALL_METRICS (MD) - textContent doesn't keep subscript
      `FRMSE`, // ALL_METRICS (MD) - textContent doesn't keep subscript
      // ΔRDF is visible:false (hidden from leaderboards, redundant with ΔvDOS/ΔADF)
      `ΔADF`, // ALL_METRICS (MD)
      `ΔvDOS`, // ALL_METRICS (MD)
      `PMAE`, // ALL_METRICS (MD) - textContent doesn't keep subscript
      `PW1`, // ALL_METRICS (MD) - textContent doesn't keep subscript
      `ΔP`, // ALL_METRICS (MD)
      `CMDS`, // ALL_METRICS (MD)
      `Speed`, // ALL_METRICS (MD)
      `Slowdown`, // ALL_METRICS (MD)
      `CDS`, // DIATOMICS_METRICS
      `Speed`, // DIATOMICS_METRICS
      `Slowdown`, // DIATOMICS_METRICS
      `E flips`, // DIATOMICS_METRICS
      `E jump`, // DIATOMICS_METRICS
      `F TV`, // DIATOMICS_METRICS
      `F flips`, // DIATOMICS_METRICS
      `F jump`, // DIATOMICS_METRICS
      `PBE ΔDe`, // DIATOMICS_METRICS
      `PBE Δr wall`, // DIATOMICS_METRICS
      `PBE Δre`, // DIATOMICS_METRICS
      `PBE Δω`, // DIATOMICS_METRICS
      `PBE E MAE`, // DIATOMICS_METRICS
      `PBE F MAE`, // DIATOMICS_METRICS
      `τ`, // DIATOMICS_METRICS
      `CPS`, // Added in assemble_row_data
    ]

    // the structural rank column (#, from show_row_numbers) is excluded by
    // header_cells() and covered by its own test below
    const header_elements = header_cells()
    const actual_core_columns = header_elements.map(header_name)

    // The default visible columns should stay intentionally curated: new default
    // columns must be added to expected_core_columns explicitly. Sorted comparison
    // ignores order but checks exact multiset (incl. duplicate Speed/Slowdown labels).
    const compare_labels = (
      label_a: string | undefined,
      label_b: string | undefined,
    ): number => (label_a ?? ``).localeCompare(label_b ?? ``)
    expect(actual_core_columns.toSorted(compare_labels)).toEqual(
      expected_core_columns.toSorted(compare_labels),
    )

    // Every column exposes its rich description through a hover/focus trigger.
    for (const header of header_elements) {
      expect(
        header.querySelector(`button[aria-haspopup]`),
        `Header ${header_name(header)} has no tooltip trigger`,
      ).not.toBeNull()
    }
  })

  it(`shows rank numbers 1..N in row order`, () => {
    mount_table()

    expect(doc_query(`thead th.row-num-col`).textContent?.trim()).toBe(`#`)
    const rank_texts = [...document.querySelectorAll(`tbody td.row-num-col`)].map((td) =>
      td.textContent?.trim(),
    )
    const n_rows = document.querySelectorAll(`tbody tr`).length
    expect(n_rows).toBeGreaterThan(0)
    expect(rank_texts).toEqual(
      Array.from({ length: n_rows }, (_, idx) => String(idx + 1)),
    )
  })

  describe(`Double-click selection functionality`, () => {
    const get_rows = () => document.querySelectorAll(`tbody tr`)
    const get_toggle = () =>
      document.querySelector<HTMLInputElement>(
        `input[aria-label="Toggle between showing only selected models and all models"]`,
      )
    const get_toggle_label = () => get_toggle()?.closest(`label`)
    const double_click_row = (row: Element) => {
      row.dispatchEvent(new MouseEvent(`dblclick`, { bubbles: true }))
    }
    const row_key = (row: Element) =>
      doc_query<HTMLAnchorElement>(`a[href^="/models/"]`, row).pathname.replace(
        `/models/`,
        ``,
      )

    it(
      `selects and deselects models on double-click with proper state management`,
      { timeout: 30_000 },
      async () => {
        mount_table({ col_filter: () => true })
        await tick() // Wait for initial render

        const initial_row = get_rows()[0]
        expect(get_rows().length).toBeGreaterThanOrEqual(2)

        // Initially no selection
        expect(get_rows()[0].classList.contains(`highlight`)).toBe(false)
        expect(get_rows()[1].classList.contains(`highlight`)).toBe(false)

        // Select first row
        double_click_row(get_rows()[0])
        await tick()
        expect(get_rows()[0].classList.contains(`highlight`)).toBe(true)
        const first_key = row_key(get_rows()[0])
        expect([...comparison.keys]).toEqual([first_key])

        // Select second row
        double_click_row(get_rows()[1])
        await tick()
        expect(get_rows()[0].classList.contains(`highlight`)).toBe(true)
        expect(get_rows()[1].classList.contains(`highlight`)).toBe(true)
        const second_key = row_key(get_rows()[1])
        expect([...comparison.keys]).toEqual([first_key, second_key])

        // Deselect first row
        double_click_row(get_rows()[0])
        await tick()
        expect(get_rows()[0].classList.contains(`highlight`)).toBe(false)
        expect(get_rows()[1].classList.contains(`highlight`)).toBe(true)
        expect([...comparison.keys]).toEqual([second_key])
        expect(get_rows()[0]).toBe(initial_row)
      },
    )

    it(
      `manages toggle visibility and count dynamically`,
      { timeout: 30_000 },
      async () => {
        mount_table({ col_filter: () => true })

        // Initially no toggle, compare button shows no count
        expect(get_toggle()).toBeNull()
        const compare_btn = doc_query<HTMLButtonElement>(`button.compare`)
        expect(compare_btn.textContent?.trim()).toBe(`Compare`)

        // Select one model
        double_click_row(get_rows()[0])
        await tick()
        expect(get_toggle()).not.toBeNull()
        expect(get_toggle_label()?.textContent).toContain(`1 selected`)
        expect(compare_btn.textContent).toContain(`Compare (1)`)

        // Select second model
        double_click_row(get_rows()[1])
        await tick()
        expect(get_toggle_label()?.textContent).toContain(`2 selected`)

        // Deselect one model
        double_click_row(get_rows()[0])
        await tick()
        expect(get_toggle_label()?.textContent).toContain(`1 selected`)

        // Deselect the remaining selected row (row 1 is still highlighted)
        double_click_row(get_rows()[1])
        await tick()

        expect(get_toggle()).toBeNull()
        expect(compare_btn.textContent?.trim()).toBe(`Compare`)

        // with nothing selected, the compare button seeds the dialog with the top 3 rows
        // in the table's current sort order (rather than opening an empty picker)
        expect(comparison.open).toBe(false)
        compare_btn.click()
        expect(comparison.open).toBe(true)
        const top_3 = () => [...get_rows()].slice(0, 3).map(row_key)
        const cps_top_3 = top_3()
        expect([...comparison.keys]).toEqual(cps_top_3)
        comparison.open = false

        // the seed follows the table's sort: alphabetical by Model gives a different top 3
        comparison.set([])
        header_cells()[0].click()
        await tick()
        compare_btn.click()
        expect(top_3()).not.toEqual(cps_top_3)
        expect([...comparison.keys]).toEqual(top_3())
        comparison.open = false

        // an existing selection is left alone
        comparison.set([cps_top_3[2]])
        compare_btn.click()
        expect([...comparison.keys]).toEqual([cps_top_3[2]])
        comparison.open = false
      },
    )

    it(`toggles models via the row context menu`, async () => {
      mount_table({ col_filter: () => true })
      await tick()
      // the toolbar with the Compare button is opted out of matterviz's hover-reveal
      expect(document.querySelector(`.table-container.leaderboard`)).not.toBeNull()

      const [row_1, row_2] = get_rows()
      // model cells hold only the link: no per-row buttons that could flash or take space
      expect(row_1.querySelector(`td[data-col="Model"] button`)).toBeNull()
      double_click_row(row_1)
      await tick()
      expect([...comparison.keys]).toEqual([row_key(row_1)])

      // right-clicking a non-link cell opens the row menu (deferred a tick past the
      // contextmenu event); links keep the browser's own menu
      const menu_items = () =>
        [...document.querySelectorAll(`menu [role="menuitem"]`)].map((el) =>
          el.textContent?.trim(),
        )
      const right_click = (el: Element) =>
        el.dispatchEvent(
          new MouseEvent(`contextmenu`, {
            bubbles: true,
            cancelable: true,
            clientX: 5,
            clientY: 5,
          }),
        )
      const flush_timers = () =>
        new Promise((resolve) => {
          setTimeout(resolve, 0)
        })
      right_click(doc_query(`td[data-col="Model"] a`, row_2))
      await flush_timers()
      await tick()
      expect(menu_items()).toEqual([])
      const cells = row_2.querySelectorAll(`td`)
      const link_event = new MouseEvent(`contextmenu`, {
        bubbles: true,
        cancelable: true,
      })
      cells[cells.length - 1].dispatchEvent(link_event)
      expect(link_event.defaultPrevented).toBe(true)
      await flush_timers()
      await tick()
      const items = menu_items()
      expect(items[0]).toContain(`to comparison`)
      expect(items[1]).toContain(`with…`)
      expect(items[2]).toContain(`model page`)
      const menu_item = (idx: number) =>
        document.querySelectorAll<HTMLButtonElement>(`menu [role="menuitem"]`)[idx]
      menu_item(0).click()
      await tick()
      expect([...comparison.keys]).toEqual([row_key(row_1), row_key(row_2)])
      expect(document.querySelector(`menu`)).toBeNull()
      // the model-page item navigates in-app
      cells[cells.length - 1].dispatchEvent(
        new MouseEvent(`contextmenu`, { bubbles: true }),
      )
      await flush_timers()
      await tick()
      menu_item(2).click()
      expect(vi.mocked(goto)).toHaveBeenCalledWith(`/models/${row_key(row_2)}`)
    })

    it(
      `filters selected rows and updates toggle labels and highlighting`,
      { timeout: 30_000 },
      async () => {
        const filters = make_table_filters()
        mount_table({ col_filter: () => true, filters })
        const initial_count = get_rows().length
        expect(initial_count).toBeGreaterThan(1)

        double_click_row(get_rows()[0])
        await tick()
        const toggle = get_toggle()
        const label = get_toggle_label()
        if (!toggle) throw new Error(`Toggle not found`)
        expect(toggle.checked).toBe(false)
        expect(label?.textContent).toContain(`Show only 1 selected`)

        toggle.click()
        await tick()
        expect(toggle.checked).toBe(true)
        expect(label?.textContent).toContain(`Show all`)
        expect(get_rows()).toHaveLength(1)
        expect(get_rows()[0].classList.contains(`highlight`)).toBe(false)
        expect(filters.url_entries).toContainEqual([`selected_only`, `1`])

        filters.read(new URLSearchParams())
        await tick()
        expect(toggle.checked).toBe(false)
        expect(get_rows()).toHaveLength(initial_count)
        filters.read(new URLSearchParams(`selected_only=1`))
        await tick()
        expect(toggle.checked).toBe(true)
        expect(get_rows()).toHaveLength(1)

        toggle.click()
        await tick()
        expect(toggle.checked).toBe(false)
        expect(label?.textContent).toContain(`Show only 1 selected`)
        expect(get_rows()).toHaveLength(initial_count)
        expect(get_rows()[0].classList.contains(`highlight`)).toBe(true)
        expect(filters.url_entries).toContainEqual([`selected_only`, ``])
      },
    )
  })

  it.each([`full_test_set`, `unique_prototypes`] as const)(
    `exports the visible ordered %s cohort directly from table data`,
    async (discovery_set) => {
      const filters = make_table_filters()
      const models = ACTIVE_MODELS.filter(
        (model) =>
          filters.matches(model) && model.metrics?.discovery?.[discovery_set]?.F1 != null,
      ).slice(0, 3)
      const keys = new Set(models.map(({ model_key }) => model_key))
      let exported_blob: Blob | undefined
      vi.spyOn(URL, `createObjectURL`).mockImplementation((blob) => {
        exported_blob = blob as Blob
        return `blob:table-export`
      })
      vi.spyOn(HTMLAnchorElement.prototype, `click`).mockImplementation(() => {})
      const copy = vi.spyOn(navigator.clipboard, `writeText`).mockResolvedValue()
      await mount_with_url(MetricsTable, `http://localhost/?sort=F1&dir=asc`, {
        props: {
          discovery_set,
          filters,
          model_filter: (model: ModelData) => keys.has(model.model_key),
          col_filter: (col: Label) => [`Model`, `F1`, `DAF`].includes(col.key),
          column_order: [`Model`, `DAF`, `F1`],
          default_sort: { column: `DAF`, dir: `desc` },
        },
      })
      await tick()
      const export_buttons = document.querySelectorAll<HTMLButtonElement>(
        `.control-buttons button[aria-haspopup="menu"]`,
      )
      expect(export_buttons).toHaveLength(1)
      expect(export_buttons[0].textContent?.trim()).toBe(`Export`)
      expect(
        doc_query(`.control-buttons`).lastElementChild?.contains(export_buttons[0]),
      ).toBe(true)
      expect(
        document.querySelector(
          `button[title^="Download visible results"], [aria-label="Copy link to this view"]`,
        ),
      ).toBeNull()
      const choose_export = async (label: string): Promise<void> => {
        export_buttons[0].click()
        await tick()
        const menu = doc_query(`menu[aria-labelledby]`)
        const items = [...menu.querySelectorAll<HTMLButtonElement>(`[role="menuitem"]`)]
        expect(items.map((item) => item.textContent?.trim())).toEqual([
          `CSV`,
          `JSON with provenance`,
          `Copy table`,
        ])
        const item = items.find((candidate) => candidate.textContent?.trim() === label)
        if (!item) throw new Error(`Missing export option: ${label}`)
        item.click()
        await tick()
        expect(document.querySelector(`menu[aria-labelledby]`)).toBeNull()
      }
      await choose_export(`CSV`)
      expect(exported_blob).toBeDefined()
      expect(exported_blob?.type).toBe(`text/csv`)
      const csv = await exported_blob?.text()
      const sorted = models.toSorted(
        (left, right) =>
          (left.metrics?.discovery?.[discovery_set]?.F1 ?? 0) -
          (right.metrics?.discovery?.[discovery_set]?.F1 ?? 0),
      )
      expect(csv?.split(`\n`)).toEqual([
        `Model,DAF,F1`,
        ...sorted.map((model) => {
          const metrics = model.metrics?.discovery?.[discovery_set]
          const excluded =
            Object.keys(model.metrics?.diatomics?.excluded_formula_reasons ?? {}).length >
            0
          return `${model.model_name}${excluded ? `*` : ``},${metrics?.DAF},${metrics?.F1}`
        }),
      ])
      await choose_export(`Copy table`)
      expect(copy).toHaveBeenLastCalledWith(csv?.replaceAll(`,`, `\t`))

      // Changing this view uses shallow URL writes: page.url intentionally stays stale.
      for (const name of [`Params`, `Training Set`, `Org`]) {
        column_checkbox(name).click()
        await tick()
      }
      const shift_click = async (id: string): Promise<void> => {
        doc_query(`th[data-col-id="${id}"]`).dispatchEvent(
          new MouseEvent(`click`, { bubbles: true, shiftKey: true }),
        )
        await tick()
      }
      await shift_click(`F1`)
      await shift_click(`model_params`)
      expect(query_param(`multi_sort`)).toBe(`-F1,-model_params`)
      const shared_url = location.href
      expect(page.url.href).not.toBe(shared_url)

      // Restore the shared view after clearing its custom columns and multiple sorts.
      const visible_ids = header_cells().map((header) => header.dataset.colId)
      const ordered_keys = row_models().map((model) => model.model_key)
      await navigate(``, `popstate`)
      expect(ranking_context()).toContain(`DAF (descending)`)
      expect(query_param(`sort`)).toBeNull()
      expect(query_param(`multi_sort`)).toBeNull()
      await navigate(new URL(shared_url).search, `popstate`)
      expect(header_cells().map((header) => header.dataset.colId)).toEqual(visible_ids)
      expect(row_models().map((model) => model.model_key)).toEqual(ordered_keys)
      expect(ranking_context()).toContain(`F1 (descending), then Params (descending)`)

      await choose_export(`JSON with provenance`)
      expect(exported_blob?.type).toBe(`application/json`)
      const exported_models = sorted.toSorted(
        (left, right) =>
          (right.metrics?.discovery?.[discovery_set]?.F1 ?? 0) -
            (left.metrics?.discovery?.[discovery_set]?.F1 ?? 0) ||
          (right.model_params ?? 0) - (left.model_params ?? 0),
      )
      expect(ordered_keys).toEqual(exported_models.map((model) => model.model_key))
      expect(JSON.parse((await exported_blob?.text()) ?? `null`)).toEqual({
        url: location.href,
        benchmark_revision: { ...BENCHMARK_REVISION, development: import.meta.env.DEV },
        exported_at: expect.any(String),
        discovery_set,
        cps_discovery_set: `unique_prototypes`,
        filters: { ...filters.as_preset, selected_only: false },
        weights: score_weight_records(),
        sort: [
          { column: `F1`, ascending: false },
          { column: `model_params`, ascending: false },
        ],
        columns: visible_ids.map((id) => expect.objectContaining({ id })),
        models: exported_models.map((model) => ({
          model_key: model.model_key,
          model_version: model.model_version,
          prediction_files: get_pred_file_urls(model),
        })),
        rows: exported_models.map((model) => ({
          Model: model.model_name,
          DAF: model.metrics?.discovery?.[discovery_set]?.DAF,
          F1: model.metrics?.discovery?.[discovery_set]?.F1,
          model_params: model.model_params,
          'Training Set': {
            datasets: model.training_sets,
            materials: model.n_training_materials,
            structures: model.n_training_structures,
          },
          Org: { logos: model.org_logos, authors: model.authors },
        })),
        references: Object.fromEntries(
          [
            `wbm_summary`,
            `wbm_initial_atoms`,
            `wbm_relaxed_atoms`,
            `wbm_dft_geo_opt_symprec_1e_2`,
            `wbm_dft_geo_opt_symprec_1e_5`,
            `phonondb_pbe_103_structures`,
            `phonondb_pbe_103_kappa_no_nac`,
            `dynamat_v1_0_md_trajectories`,
            `diatomics_dft_reference`,
          ].map((key) => {
            const entry = data_files[key]
            if (typeof entry === `string`)
              throw new Error(`Invalid reference file: ${key}`)
            const { url, path, md5 } = entry
            return [key, { url, path, md5 }]
          }),
        ),
      })

      // Every format also honors selected-only filtering, including an empty selection.
      filters.show_selected_only = true
      comparison.keys.add(exported_models[0].model_key)
      await tick()
      await choose_export(`Copy table`)
      expect(copy.mock.lastCall?.[0].split(`\n`)).toHaveLength(2)
      expect(copy.mock.lastCall?.[0]).toContain(exported_models[0].model_name)
      comparison.keys.clear()
      await tick()
      await choose_export(`CSV`)
      expect((await exported_blob?.text())?.split(`\n`)).toHaveLength(1)
      await choose_export(`JSON with provenance`)
      expect(JSON.parse((await exported_blob?.text()) ?? `null`)).toMatchObject({
        filters: { selected_only: true },
        models: [],
        rows: [],
      })
    },
  )

  describe(`Column Reordering`, () => {
    it.each([
      [`?columns=Model,F1&column_order=F1,Model`, [`F1`, `Model`]],
      [`?columns=Model,F1,F1&column_order=F1,F1,Model`, [`F1`, `Model`]],
      [
        `?columns=Model,missing&column_order=missing&multi_sort=-missing`,
        [`Model`, `DAF`, `F1`],
      ],
      [`?columns=none`, []],
    ])(`restores and validates custom columns from %s`, async (query, expected) => {
      await mount_with_url(MetricsTable, `http://localhost/${query}`, {
        props: {
          model_filter: (model: ModelData) => model === ACTIVE_MODELS[0],
          col_filter: (col: Label) => [`Model`, `F1`, `DAF`].includes(col.key),
          column_order: [`Model`, `DAF`, `F1`],
        },
      })
      expect(header_names()).toEqual(expected)
      expect(query_param(`multi_sort`)).toBeNull()
      await navigate(``, `popstate`)
      expect(header_names()).toEqual([`Model`, `DAF`, `F1`])
      expect(query_param(`columns`)).toBeNull()
      expect(query_param(`column_order`)).toBeNull()
    })

    it.each([false, true])(
      `keeps labels and explicit columns consistent across mobile=%s and resizing`,
      async (initial_mobile) => {
        let is_mobile = initial_mobile
        const media = Object.assign(new EventTarget(), {
          matches: is_mobile,
          media: `(max-width: 600px)`,
          onchange: null,
          addListener: () => {},
          removeListener: () => {},
        }) as MediaQueryList
        Object.defineProperty(media, `matches`, { get: () => is_mobile })
        const match_media = window.matchMedia
        vi.spyOn(window, `matchMedia`).mockImplementation((query) =>
          query === media.media ? media : match_media(query),
        )
        const resize = async (mobile: boolean): Promise<void> => {
          is_mobile = mobile
          media.dispatchEvent(new Event(`change`))
          await tick()
        }
        const diatomics_keys = new Set([
          `Model`,
          ...Object.values(DIATOMICS_METRICS).map((col) => col.key),
          `model_params`,
        ])
        mount_table({
          model_filter: (model: ModelData) => model === ACTIVE_MODELS[0],
          col_filter: (col: Label) => diatomics_keys.has(col.key),
          default_sort: { column: `diatomics_combined_score`, dir: `desc` },
        })
        // The first client render must match server HTML before onMount compacts it.
        expect(header_cells().length).toBeGreaterThan(5)
        await tick()
        await resize(true)
        expect(header_names()).toEqual([`Model`, `τ`, `E flips`, `CDS`, `Params`])
        expect(header_cells().map((header) => header.dataset.colId)).toEqual(
          query_param(`columns`)?.split(`,`),
        )
        const mobile_query = location.search
        await resize(false)
        expect(header_cells().length).toBeGreaterThan(5)
        expect(query_param(`columns`)).toBeNull()

        // The shared mobile view also restores those exact columns on desktop.
        await navigate(mobile_query, `popstate`)
        expect(header_names()).toEqual([`Model`, `τ`, `E flips`, `CDS`, `Params`])
        column_checkbox(`Org`).click()
        await tick()
        const custom_headers = header_names()
        expect(custom_headers).toContain(`Org`)
        await resize(true)
        expect(header_names()).toEqual(custom_headers)
        await resize(false)
        expect(header_names()).toEqual(custom_headers)
      },
    )

    it(`initializes all columns and displays visible columns in column_order`, async () => {
      const state = { column_order: [] as string[] }
      mount_table({
        get column_order() {
          return state.column_order
        },
        set column_order(val) {
          state.column_order = val
        },
        col_filter: (col: Label) => [`Model`, `F1`, `DAF`].includes(col.key ?? col.label),
      })
      await tick()

      // After mounting, column_order should be initialized with ALL columns
      // (not just visible ones - the visible filter is separate)
      expect(state.column_order.length).toBeGreaterThan(10)
      expect(state.column_order).toContain(`Model`)
      expect(state.column_order).toContain(`F1`)
      expect(state.column_order).toContain(`DAF`)

      expect(header_names()).toStrictEqual([`Model`, `F1`, `DAF`])

      const f1_idx = state.column_order.indexOf(`F1`)
      const daf_idx = state.column_order.indexOf(`DAF`)
      expect(f1_idx).toBeGreaterThanOrEqual(0)
      expect(daf_idx).toBeGreaterThanOrEqual(0)

      const headers = header_names()
      expect(headers[0]).toBe(`Model`)

      // F1 and DAF should appear in the order specified by column_order
      const visible_f1_pos = headers.indexOf(`F1`)
      const visible_daf_pos = headers.indexOf(`DAF`)
      expect(
        f1_idx < daf_idx
          ? visible_f1_pos < visible_daf_pos
          : visible_f1_pos > visible_daf_pos,
      ).toBe(true)
    })

    it.each([
      { columns: [`Model`, `F1`, `DAF`], name: `basic columns` },
      { columns: [`Model`, `F1`, `DAF`, `CPS`], name: `with CPS` },
    ])(`maintains Model column first with $name`, async ({ columns }) => {
      mount_table({
        col_filter: (col: Label) => columns.includes(col.key ?? col.label),
      })
      await tick()

      const headers = header_cells()
      expect(headers[0].textContent?.split(` `)[0]).toBe(`Model`)
      expect(headers[0].classList.contains(`sticky-col`)).toBe(true)
    })

    it(`preserves column_order when toggling column visibility`, async () => {
      const state = {
        col_filter: (col: Label) =>
          [`Model`, `F1`, `DAF`, `CPS`].includes(col.key ?? col.label),
        column_order: [] as string[],
      }

      mount_table({
        get col_filter() {
          return state.col_filter
        },
        get column_order() {
          return state.column_order
        },
        set column_order(val) {
          state.column_order = val
        },
      })
      await tick()

      const initial_order = [...state.column_order]
      expect(initial_order.length).toBeGreaterThan(10)
      const [f1_idx, daf_idx, cps_idx] = [
        initial_order.indexOf(`F1`),
        initial_order.indexOf(`DAF`),
        initial_order.indexOf(`CPS`),
      ]

      state.col_filter = (col: Label) =>
        [`Model`, `F1`, `DAF`].includes(col.key ?? col.label)
      await tick()

      expect(state.column_order).toHaveLength(initial_order.length)
      expect(state.column_order).toContain(`CPS`) // Still in order, just not visible

      // Positions should be unchanged
      expect(state.column_order.indexOf(`F1`)).toBe(f1_idx)
      expect(state.column_order.indexOf(`DAF`)).toBe(daf_idx)
      expect(state.column_order.indexOf(`CPS`)).toBe(cps_idx)
    })

    it(`sets columns as draggable without initial drag state`, () => {
      mount_table({
        col_filter: (col: Label) => [`Model`, `F1`, `DAF`].includes(col.key ?? col.label),
      })

      // header_cells() excludes the rank (#) column, which is structural and
      // deliberately not draggable
      const headers = header_cells()
      expect(headers.length).toBeGreaterThan(0)
      headers.forEach((header) => {
        expect(header.getAttribute(`draggable`)).toBe(`true`)
        expect(header.getAttribute(`aria-dropeffect`)).toBe(`move`)
        // no drag state classes before any drag interaction
        expect(header.classList.contains(`dragging`)).toBe(false)
        expect(header.classList.contains(`drag-over`)).toBe(false)
      })
    })
  })
})
