import { HYPERPARAMS } from '$lib/labels'
import { ACTIVE_MODELS, make_table_filters } from '$lib/models.svelte'
import MetricsTable from '$lib/table/MetricsTable.svelte'
import type { Label, ModelData } from '$lib/types'
import { tick } from 'svelte'
import { describe, expect, it } from 'vitest'
import { doc_query, mount } from '../index'

// all header cells except the structural rank (#) column HeatmapTable renders
// for show_row_numbers
const header_cells = () => [
  ...document.querySelectorAll<HTMLTableCellElement>(`th:not(.row-num-col)`),
]
const header_name = (header: HTMLTableCellElement) =>
  header.querySelector(`.header-label`)?.textContent?.trim()
const header_names = () => header_cells().map(header_name)

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
  const parse_integer_sort_value = (cell: Element): number | null => {
    const sort_value = cell.getAttribute(`data-sort-value`)
    if (!sort_value) return null

    const parsed_value = Math.trunc(Number(sort_value))
    return Number.isNaN(parsed_value) ? null : parsed_value
  }

  it(`renders with default props`, async () => {
    mount(MetricsTable, {
      target: document.body,
      props: { col_filter: () => true },
    })

    // Check table structure
    const table = document.querySelector(`table`)
    expect(table).toBeDefined()
    expect(table?.querySelector(`thead`)).toBeDefined()
    expect(table?.querySelector(`tbody`)).toBeDefined()
    const table_container = doc_query(`.table-container`)
    expect(
      table_container.style.getPropertyValue(`--heatmap-sticky-cell-odd-bg`),
    ).toContain(`linear-gradient`)
    expect(table_container.style.getPropertyValue(`--heatmap-row-num-padding-left`)).toBe(
      `0`,
    )

    // Check essential columns are present (with sort indicators)
    const header_texts = header_cells().map((h) => h.textContent?.trim())
    const required_cols = [
      `Model`,
      `CPS ↑`, // active sort column has indicator
      `F1`,
      `DAF`,
      `Training Set`,
      `Params`,
      `Targets`,
      `Links`,
    ]

    // Make sure each required column is present
    for (const col of required_cols) {
      expect(header_texts).toContain(col)
    }

    // Model stays first and Org is a regular metadata column at the far right.
    expect(header_texts[0]).toBe(`Model`)
    expect(header_texts.at(-1)).toBe(`Org`)
    const metric_order = [`CPS ↑`, `F1`, `DAF`].map((col) => header_texts.indexOf(col))
    expect(metric_order).toStrictEqual([...metric_order].toSorted((n1, n2) => n1 - n2))

    // Test prediction files dropdown interaction
    const pred_files_button = doc_query<HTMLButtonElement>(`tbody .pred-files-btn`)
    expect(pred_files_button).toBeDefined() // Ensure at least one button exists

    // Dropdown should not exist initially
    expect(document.querySelector(`.pred-files-dropdown`)).toBeNull()

    // Click the button
    pred_files_button.click()
    await tick() // Wait for state update and render

    // Dropdown should now exist
    let dropdown = document.querySelector(`.pred-files-dropdown`)
    expect(dropdown).toBeDefined()
    expect(dropdown?.textContent).toContain(`Files for`)

    // Simulate clicking outside (assuming click_outside works correctly)
    // We can simulate this by clicking the body or another element
    document.body.click()
    await tick()

    // Dropdown should be gone after clicking outside
    dropdown = document.querySelector(`.pred-files-dropdown`)
    expect(dropdown).toBeNull()

    // Test closing with Escape key
    pred_files_button.click() // Reopen dropdown
    await tick()
    dropdown = document.querySelector(`.pred-files-dropdown`)
    expect(dropdown).toBeDefined()

    // Dispatch Escape keydown event
    globalThis.dispatchEvent(new KeyboardEvent(`keydown`, { key: `Escape` }))
    await tick()

    // Dropdown should be gone
    dropdown = document.querySelector(`.pred-files-dropdown`)
    expect(dropdown).toBeNull()
  })

  it(`renders Org as a regular rightmost metadata column`, async () => {
    mount(MetricsTable, {
      target: document.body,
      props: { col_filter: () => true },
    })
    await tick()

    const org_cell = doc_query(`td[data-col="Org"]`)
    const org_preview = doc_query(`td[data-col="Org"] .org-preview`)
    const headers = header_cells()
    const org_header = headers.at(-1)
    if (!org_header) throw new Error(`Org column header not found`)

    expect(org_header?.textContent?.trim()).toBe(`Org`)
    expect(org_header.getAttribute(`title`)).toBeNull()
    expect(org_header.querySelector(`.header-label`)).not.toBeNull()
    expect(org_preview.classList.contains(`org-preview`)).toBe(true)
    expect(org_cell.getAttribute(`style`)).not.toContain(`min-width:`)
  })

  it(`renders header tooltips on inner labels`, async () => {
    mount(MetricsTable, {
      target: document.body,
      props: { col_filter: () => true },
    })
    await tick()

    const cps_header = header_cells().find((header) =>
      header.textContent?.trim().startsWith(`CPS`),
    )
    if (!cps_header) throw new Error(`CPS column header not found`)

    cps_header.dispatchEvent(new MouseEvent(`mouseover`, { bubbles: true }))
    await tick()
    expect(cps_header.getAttribute(`title`)).toBeNull()
    expect(cps_header.querySelector(`.header-label`)).not.toBeNull()
  })

  it(`toggles metadata columns`, () => {
    // Keys used by col_filter (col.key ?? col.label)
    const metadata_keys = new Set([`Training Set`, `Targets`, `benchmark_added`, `Links`])
    // Labels displayed in table headers
    const metadata_labels = [`Training Set`, `Targets`, `Date Added`, `Links`]
    const col_filter = (_col: Label) => true // show all columns initially
    mount(MetricsTable, { target: document.body, props: { col_filter } })
    // Check metadata columns are visible initially
    let header_texts = header_names()
    expect(header_texts).toStrictEqual(expect.arrayContaining(metadata_labels))

    // Create a new instance that hides metadata columns
    document.body.innerHTML = ``
    mount(MetricsTable, {
      target: document.body,
      props: {
        col_filter: (col: Label) => !metadata_keys.has(col.key ?? col.label),
      },
    })

    // Check metadata columns are hidden
    header_texts = header_names()

    // Each metadata column label should be hidden
    for (const col of metadata_labels) {
      expect(header_texts).not.toContain(col)
    }

    // Check metric columns are still visible (strip sort indicators)
    const metric_cols = [`CPS`, `F1`, `DAF`, `Prec`, `Acc`]
    for (const col of metric_cols) {
      expect(header_texts).toContain(col)
    }
  })

  it(`filters specified columns`, () => {
    const col_filter = (col: Label) => ![`F1`, `DAF`].includes(col.key ?? col.label)
    mount(MetricsTable, {
      target: document.body,
      props: { col_filter },
    })

    const header_texts = header_names()

    // Check hidden columns
    expect(header_texts).not.toContain(`F1`)
    expect(header_texts).not.toContain(`DAF`)

    // Check other columns still visible
    expect(header_texts).toContain(`Model`)
    expect(header_texts).toContain(`CPS`)
    expect(header_texts).toContain(`Prec`)
    expect(header_texts).toContain(`Acc`)
  })

  it(
    `hides energy-only models by default via the targets filter`,
    {
      timeout: 30_000,
    },
    async () => {
      // default filters require force prediction, hiding energy-only models
      mount(MetricsTable, { target: document.body })
      await tick()
      const rows_without_energy = document.querySelectorAll(`tbody tr`).length

      // clearing the targets filter shows them
      document.body.innerHTML = ``
      const show_all = all_targets_filters()
      mount(MetricsTable, { target: document.body, props: { filters: show_all } })
      await tick()
      const rows_with_energy = document.querySelectorAll(`tbody tr`).length

      expect(rows_without_energy).toBe(visible_row_count())
      expect(rows_with_energy).toBe(visible_row_count(show_all.matches))
      expect(rows_with_energy).toBeGreaterThan(rows_without_energy)
    },
  )

  it(`filters models based on model_filter prop`, () => {
    // First test: show no models
    const no_model_filter = (_model: ModelData) => false
    mount(MetricsTable, {
      target: document.body,
      props: {
        model_filter: no_model_filter,
      },
    })

    const data_rows = document.querySelectorAll(`tbody tr`)
    // HeatmapTable may render a "no data" placeholder row when empty
    expect(data_rows.length).toBeLessThanOrEqual(1)

    // Second test: show all models
    document.body.innerHTML = ``
    mount(MetricsTable, {
      target: document.body,
      props: {
        model_filter: () => true,
      },
    })

    const all_rows = document.querySelectorAll(`tbody tr`).length
    expect(all_rows).toBe(visible_row_count())

    // Third test: show specific models (e.g., only models with CHG in name)
    document.body.innerHTML = ``
    mount(MetricsTable, {
      target: document.body,
      props: {
        model_filter: (model: ModelData) => model.model_name.includes(`CHG`),
      },
    })

    const filtered_rows = document.querySelectorAll(`tbody tr`)
    const default_matches = make_table_filters().matches
    expect(filtered_rows).toHaveLength(
      visible_row_count(
        (model) => default_matches(model) && model.model_name.includes(`CHG`),
      ),
    )
    expect(filtered_rows.length).toBeLessThan(all_rows)

    // Verify that filtered rows actually contain CHG
    filtered_rows.forEach((row) => {
      const model_cell = row.querySelector(`td[data-col="Model"]`)
      expect(model_cell?.textContent).toContain(`CHG`)
    })
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
    mount(MetricsTable, { target: document.body, props: { col_filter } })

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
      mount(MetricsTable, {
        target: document.body,
        props: { model_filter, col_filter },
      })

      expect(header_names()).toStrictEqual(expected_headers)

      const rows = document.querySelectorAll(`tbody tr`)
      rows.forEach((row) => {
        const model_cell = row.querySelector(`td[data-col="Model"]`)
        expect(model_cell?.textContent).toContain(expected_model_match)
      })
    },
  )

  it(`marks models with excluded metric samples`, async () => {
    document.body.innerHTML = ``
    mount(MetricsTable, {
      target: document.body,
      props: {
        model_filter: (model: ModelData) => model.model_name === `AlphaNet-v1-OAM`,
        col_filter: (col: Label) => col.label === `Model`,
      },
    })
    await tick()

    const model_cell = doc_query(`td[data-col="Model"]`)
    const marker = model_cell.querySelector(
      `span[data-original-title="Diatomics metrics exclude He-He due to exploding errors"]`,
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
        data_col: `Date Added`,
      },
      { col_key: `Training Set`, header: `Training Set`, data_col: `Training Set` },
      {
        col_key: HYPERPARAMS.model_params.key,
        header: `Params`,
        data_col: HYPERPARAMS.model_params.label,
      },
    ])(
      `sorts $header numerically via data-sort-value`,
      async ({ col_key, header, data_col }) => {
        mount(MetricsTable, {
          target: document.body,
          props: {
            col_filter: (col: Label) => [`Model`, col_key].includes(col.key ?? col.label),
          },
        })

        const sort_header = header_cells().find((th) => th.textContent?.includes(header))
        if (!sort_header) throw new Error(`${header} column not found`)

        const cell_values = () =>
          [...document.querySelectorAll(`td[data-col="${data_col}"]`)]
            .map((cell) => {
              const sortable = cell.querySelector(`[data-sort-value]`)
              if (sortable) return parse_integer_sort_value(sortable)
              const timestamp = Date.parse(cell.textContent?.trim() ?? ``)
              return Number.isNaN(timestamp) ? null : timestamp
            })
            .filter((val) => val !== null)

        sort_header.click()
        await tick()

        const values = cell_values()
        expect(values.length).toBeGreaterThan(1)
        // rows must be fully sorted (either direction) after the click
        const ascending = [...values].toSorted((val_1, val_2) => val_1 - val_2)
        expect([ascending, ascending.toReversed()]).toContainEqual(values)

        // second click toggles sort direction
        sort_header.click()
        await tick()
        expect(cell_values()).toStrictEqual(values.toReversed())
      },
    )

    it(`properly handles HTML content in cells without using it for data-sort-value`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: {
          col_filter: (col: Label) => [`Model`, `Training Set`].includes(col.label),
        },
      })

      await tick()

      // Find all Training Set cells
      const training_set_cells = [
        ...document.querySelectorAll(`td[data-col="Training Set"]`),
      ]

      // Find cells with HTML content (looking for cells containing spans with tooltips)
      const html_cells = training_set_cells.filter(
        (cell) => cell.innerHTML.includes(`<span`) && cell.innerHTML.includes(`title=`),
      )

      // Ensure we found some cells with HTML content
      expect(html_cells.length).toBeGreaterThan(0)

      // Check that data-sort-value attribute on the td is not the full HTML
      html_cells.forEach((cell) => {
        const data_sort_value = cell.getAttribute(`data-sort-value`)

        // The data-sort-value should not contain HTML tags if present
        expect(
          data_sort_value == null ||
            (!/[<>]/.test(data_sort_value) && !data_sort_value.includes(`span`)),
        ).toBe(true)

        // The inner span should have its own data-sort-value
        const inner_span = cell.querySelector(`span[data-sort-value]`)
        const span_sort_value = inner_span?.getAttribute(`data-sort-value`)
        expect(
          inner_span === null ||
            (span_sort_value !== null && !Number.isNaN(Number(span_sort_value))),
        ).toBe(true)
      })

      // Verify tooltips are preserved on spans within cells
      const cells_with_tooltips = training_set_cells.filter(
        (cell) => cell.querySelector(`span[data-original-title]`) !== null,
      )

      expect(cells_with_tooltips.length).toBeGreaterThan(0)
      cells_with_tooltips.forEach((cell) => {
        const span = cell.querySelector(`span[data-original-title]`)
        expect(span?.getAttribute(`data-original-title`)).not.toBeNull()
      })
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
        mount(MetricsTable, { target: document.body, props })

        // Find Model column header
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

        model_header.click() // Click to sort (ascending A-Z)
        await tick()

        // Get model names after first sort
        const sorted_model_names = get_model_names()

        expect(sorted_model_names).toHaveLength(
          visible_row_count(
            `filters` in props && props.filters ? props.filters.matches : undefined,
          ),
        )

        // Verify sorted in some alphabetical order (ascending or descending)
        const ascending = [...sorted_model_names].toSorted((a, b) => a.localeCompare(b))
        const is_ascending =
          JSON.stringify(sorted_model_names) === JSON.stringify(ascending)
        const is_descending =
          JSON.stringify(sorted_model_names) === JSON.stringify(ascending.toReversed())
        expect(is_ascending || is_descending).toBe(true)

        // Click again to reverse sort direction
        model_header.click()
        await tick()

        const reverse_sorted_model_names = get_model_names()
        // Second click should reverse the sort direction
        expect(reverse_sorted_model_names).toStrictEqual(
          is_ascending ? ascending.toReversed() : ascending,
        )
      },
    )

    it(`prevents sorting of unsortable Links column`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: {
          col_filter: (col: Label) =>
            [`Model`, `CPS`, `Links`].includes(col.key ?? col.label),
        },
      })

      // Find CPS and Links column headers
      const headers = header_cells()
      const cps_header = headers.find((h) => h.textContent?.includes(`CPS`))
      if (!cps_header) throw new Error(`CPS column not found`)
      const links_header = headers.find((h) => h.textContent?.includes(`Links`))
      if (!links_header) throw new Error(`Links column not found`)

      // Verify Links header has not-sortable class
      expect(links_header.classList.contains(`not-sortable`)).toBe(true)

      // Sort by CPS first (to establish a known order)
      cps_header.click()
      await tick()

      // Get model names in current order
      const initial_models = [...document.querySelectorAll(`td[data-col="Model"]`)].map(
        (cell) => cell.textContent,
      )

      // Try to sort by Links
      links_header.click()
      await tick()

      // Get model names after clicking Links
      const after_links_click_models = [
        ...document.querySelectorAll(`td[data-col="Model"]`),
      ].map((cell) => cell.textContent)

      // Order should not change
      expect(after_links_click_models).toStrictEqual(initial_models)
    })
  })

  describe(`Links Column`, () => {
    it(`renders external links with proper attributes`, async () => {
      const col_filter = (col: Label) => [`Model`, `Links`].includes(col.label)
      mount(MetricsTable, {
        target: document.body,
        props: { col_filter },
      })

      await tick() // Wait for component to process data

      // Find all links cells
      const links_cells = [...document.querySelectorAll(`td[data-col="Links"]`)]
      expect(links_cells).toHaveLength(visible_row_count())

      // Check that rows have links (at least some should)
      let rows_with_links = 0
      for (const cell of links_cells) {
        const links = [...cell.querySelectorAll(`a`)]
        if (links.length > 1) rows_with_links++

        // Check each link has proper attributes
        for (const link of links) {
          expect(link.getAttribute(`target`)).toBe(`_blank`)
          expect(link.getAttribute(`rel`)).toBe(`noopener noreferrer`)

          const title = link.getAttribute(`data-original-title`)
          const href = link.getAttribute(`href`)
          expect(title).not.toBeNull()
          expect(title).not.toBe(``)
          expect(href).not.toBeNull()
          expect(href).not.toBe(``)

          // Each link should have an SVG icon - use a more general selector
          const svg = link.querySelector(`svg`)
          expect(svg).not.toBeNull()
        }
      }

      // At least half of rows should have multiple links
      expect(rows_with_links).toBeGreaterThan(links_cells.length / 2)
    })

    it(`shows icon-unavailable for missing links`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: {
          col_filter: (col: Label) => [`Model`, `Links`].includes(col.label),
        },
      })

      await tick() // Wait for component to process data

      // Find all links cells
      const links_cells = [...document.querySelectorAll(`td[data-col="Links"]`)]

      const missing_icon_titles = links_cells.flatMap((cell) =>
        [...cell.querySelectorAll(`span[data-original-title$="not available"] svg`)].map(
          (icon) => icon.closest(`span`)?.getAttribute(`data-original-title`),
        ),
      )

      // Note: this might fail if all models have all links, which is unlikely
      expect(missing_icon_titles.length).toBeGreaterThan(0)
      expect(missing_icon_titles.every((title) => title?.match(/not available/))).toBe(
        true,
      )
    })

    it(`renders prediction files button`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: {
          col_filter: (col: Label) => [`Model`, `Links`].includes(col.label),
        },
      })

      await tick() // Wait for component to process data

      // Find all pred_files buttons (every row renders one)
      const pred_file_buttons = [...document.querySelectorAll(`.pred-files-btn`)]
      expect(pred_file_buttons).toHaveLength(visible_row_count())

      // Check button attributes
      for (const button of pred_file_buttons) {
        expect(button.getAttribute(`aria-label`)).toBe(`Download model prediction files`)

        // Check for the SVG icon
        const svg = button.querySelector(`svg`)
        expect(svg).not.toBeNull()
      }
    })

    it(`renders consistent link order and specific icons`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: {
          col_filter: (col: Label) => [`Model`, `Links`].includes(col.label),
        },
      })

      await tick() // Wait for component to process data

      // Find all links cells with at least one link
      const links_cells = [...document.querySelectorAll(`td[data-col="Links"]`)].filter(
        (cell) => cell.querySelectorAll(`a`).length > 0,
      )

      // There should be at least one cell with links
      expect(links_cells.length).toBeGreaterThan(0)

      // Get the first cell with links to use as reference
      const reference_cell = links_cells[0]
      const reference_links = [...reference_cell.querySelectorAll(`a`)]

      // Check that links have SVG icons (now using Icon component)
      const reference_icons = reference_links.map((link) => {
        const svg = link.querySelector(`svg`)
        return svg !== null
      })

      // Verify that links have SVG icons
      expect(reference_icons.length).toBeGreaterThan(0)
      expect(reference_icons.some(Boolean)).toBe(true)

      // Check that the order of links is consistent across cells
      // (not checking all cells as some might have missing links)
      const second_cell = links_cells[1]
      const second_links = second_cell ? [...second_cell.querySelectorAll(`a`)] : []
      const second_icons = second_links.map((link) => {
        const svg = link.querySelector(`svg`)
        return svg !== null
      })
      expect(
        links_cells.length <= 1 ||
          second_links.length !== reference_links.length ||
          second_icons.every(Boolean),
      ).toBe(true)
    })
  })

  describe(`Heatmap Toggle Interaction`, () => {
    // happy-dom renders of the full-column table are slow in CI
    it(
      `toggles heatmap colors via TableControls checkbox`,
      { timeout: 30_000 },
      async () => {
        mount(MetricsTable, {
          target: document.body,
          props: {},
        })
        await tick() // Wait for initial render

        // Find the heatmap toggle checkbox within TableControls
        const heatmap_checkbox = document.querySelector<HTMLInputElement>(
          `input[type="checkbox"][aria-label="Toggle heatmap colors"]`,
        )

        expect(heatmap_checkbox).not.toBeNull()
        if (!heatmap_checkbox) return // Type guard

        // Initially, heatmap should be on (default)
        expect(heatmap_checkbox.checked).toBe(true)

        // Click the checkbox to turn heatmap off
        heatmap_checkbox.click()
        await tick()

        expect(heatmap_checkbox.checked).toBe(false)

        // Click again to turn heatmap back on
        heatmap_checkbox.click()
        await tick()

        expect(heatmap_checkbox.checked).toBe(true)
      },
    )
  })

  it(`renders the correct default columns`, () => {
    mount(MetricsTable, { target: document.body })

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
      // Recall is visible: false
      `TNR`, // ALL_METRICS (Discovery)
      `TPR`, // ALL_METRICS (Discovery)
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

    // Header tooltip content is attached to inner labels so HeatmapTable's
    // generic title-based tooltip doesn't flash below before our desired top placement.
    header_elements.forEach((th) => {
      const title = th.getAttribute(`title`)
      expect(title, `Header ${th.textContent} has stale title`).toBeNull()
      expect(
        th.querySelector(`.header-label`),
        `Header ${th.textContent} has no tooltip label`,
      ).not.toBeNull()
    })
  })

  it(`shows rank numbers 1..N in row order`, () => {
    mount(MetricsTable, { target: document.body })

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

    it(
      `selects and deselects models on double-click with proper state management`,
      { timeout: 30_000 },
      async () => {
        mount(MetricsTable, {
          target: document.body,
          props: { col_filter: () => true },
        })
        await tick() // Wait for initial render

        // Use fresh row references each time to avoid stale DOM after re-renders
        expect(get_rows().length).toBeGreaterThanOrEqual(2)

        // Initially no selection
        expect(get_rows()[0].classList.contains(`highlight`)).toBe(false)
        expect(get_rows()[1].classList.contains(`highlight`)).toBe(false)

        // Select first row
        double_click_row(get_rows()[0])
        await tick()
        expect(get_rows()[0].classList.contains(`highlight`)).toBe(true)

        // Select second row
        double_click_row(get_rows()[1])
        await tick()
        expect(get_rows()[0].classList.contains(`highlight`)).toBe(true)
        expect(get_rows()[1].classList.contains(`highlight`)).toBe(true)

        // Deselect first row
        double_click_row(get_rows()[0])
        await tick()
        expect(get_rows()[0].classList.contains(`highlight`)).toBe(false)
        expect(get_rows()[1].classList.contains(`highlight`)).toBe(true)
      },
    )

    it(
      `manages toggle visibility and count dynamically`,
      { timeout: 30_000 },
      async () => {
        mount(MetricsTable, {
          target: document.body,
          props: { col_filter: () => true },
        })

        // Initially no toggle
        expect(get_toggle()).toBeNull()

        // Select one model
        double_click_row(get_rows()[0])
        await tick()
        expect(get_toggle()).not.toBeNull()
        expect(get_toggle_label()?.textContent).toContain(`1 selected`)

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
      },
    )

    it(
      `toggles filter state and updates UI labels correctly`,
      { timeout: 30_000 },
      async () => {
        mount(MetricsTable, {
          target: document.body,
          props: { col_filter: () => true },
        })

        // Select a model to make toggle visible
        double_click_row(get_rows()[0])
        await tick()

        const toggle = get_toggle()
        const label = get_toggle_label()
        expect(toggle).not.toBeNull()
        if (!toggle) return // Type guard

        // Test toggle states
        expect(toggle.checked).toBe(false)
        expect(label?.textContent).toContain(`Show only 1 selected`)

        toggle.click()
        await tick()
        expect(toggle.checked).toBe(true)
        expect(label?.textContent).toContain(`Show all`)

        toggle.click()
        await tick()
        expect(toggle.checked).toBe(false)
        expect(label?.textContent).toContain(`Show only 1 selected`)
      },
    )

    it(
      `filters rows and manages styling based on filter state`,
      { timeout: 30_000 },
      async () => {
        mount(MetricsTable, {
          target: document.body,
          props: { col_filter: () => true },
        })

        const initial_count = get_rows().length
        expect(initial_count).toBeGreaterThan(1)

        // Select first row
        double_click_row(get_rows()[0])
        await tick()

        // Enable filter
        const toggle_enable = get_toggle()
        if (!toggle_enable) throw new Error(`Toggle not found`)
        toggle_enable.click()
        await tick()

        // Should show only selected row
        expect(get_rows()).toHaveLength(1)
        expect(get_rows()[0].classList.contains(`highlight`)).toBe(false) // No highlight when filtering

        // Disable filter
        const toggle_disable = get_toggle()
        if (!toggle_disable) throw new Error(`Toggle not found`)
        toggle_disable.click()
        await tick()

        // Should show all rows with highlight
        expect(get_rows()).toHaveLength(initial_count)
        expect(get_rows()[0].classList.contains(`highlight`)).toBe(true)
      },
    )
  })

  describe(`Column Reordering`, () => {
    it(`initializes column_order with all columns, not just visible ones`, async () => {
      const state = { column_order: [] as string[] }
      mount(MetricsTable, {
        target: document.body,
        props: {
          get column_order() {
            return state.column_order
          },
          set column_order(val) {
            state.column_order = val
          },
          col_filter: (col: Label) =>
            [`Model`, `F1`, `DAF`].includes(col.key ?? col.label),
        },
      })
      await tick()

      // After mounting, column_order should be initialized with ALL columns
      // (not just visible ones - the visible filter is separate)
      expect(state.column_order.length).toBeGreaterThan(10)
      expect(state.column_order).toContain(`Model`)
      expect(state.column_order).toContain(`F1`)
      expect(state.column_order).toContain(`DAF`)

      expect(header_names()).toStrictEqual([`Model`, `F1`, `DAF`])
    })

    it.each([
      { columns: [`Model`, `F1`, `DAF`], name: `basic columns` },
      { columns: [`Model`, `F1`, `DAF`, `CPS`], name: `with CPS` },
    ])(`maintains Model column first with $name`, async ({ columns }) => {
      mount(MetricsTable, {
        target: document.body,
        props: {
          col_filter: (col: Label) => columns.includes(col.key ?? col.label),
        },
      })
      await tick()

      const headers = header_cells()
      expect(headers[0].textContent?.split(` `)[0]).toBe(`Model`)
      expect(headers[0].classList.contains(`sticky-col`)).toBe(true)
    })

    it(`respects column_order for visible column display order`, async () => {
      const state = { column_order: [] as string[] }
      mount(MetricsTable, {
        target: document.body,
        props: {
          get column_order() {
            return state.column_order
          },
          set column_order(val) {
            state.column_order = val
          },
          col_filter: (col: Label) =>
            [`Model`, `F1`, `DAF`].includes(col.key ?? col.label),
        },
      })
      await tick()

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

    it(`preserves column_order when toggling column visibility`, async () => {
      const state = {
        col_filter: (col: Label) =>
          [`Model`, `F1`, `DAF`, `CPS`].includes(col.key ?? col.label),
        column_order: [] as string[],
      }

      mount(MetricsTable, {
        target: document.body,
        props: {
          get col_filter() {
            return state.col_filter
          },
          get column_order() {
            return state.column_order
          },
          set column_order(val) {
            state.column_order = val
          },
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
      mount(MetricsTable, {
        target: document.body,
        props: {
          col_filter: (col: Label) =>
            [`Model`, `F1`, `DAF`].includes(col.key ?? col.label),
        },
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
