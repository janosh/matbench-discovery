import { HYPERPARAMS } from '$lib/labels'
import { model_is_compliant, MODELS } from '$lib/models.svelte'
import MetricsTable from '$lib/table/MetricsTable.svelte'
import type { Label, ModelData } from '$lib/types'
import { tick } from 'svelte'
import { describe, expect, it } from 'vitest'
import { doc_query, mount } from '../index'

// expected table row count for the default unique_prototypes discovery set,
// mirroring MetricsTable's model filters
const visible_row_count = (
  { energy = false, non_compliant = true } = {},
  extra_filter: (model: ModelData) => boolean = () => true,
) =>
  MODELS.filter(
    (model) =>
      (energy || model.targets !== `E`) &&
      (non_compliant || model_is_compliant(model)) &&
      typeof model.metrics?.discovery === `object` &&
      model.metrics.discovery.unique_prototypes &&
      extra_filter(model),
  ).length

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
      props: { col_filter: () => true, show_non_compliant: true },
    })

    // Check table structure
    const table = document.querySelector(`table`)
    expect(table).toBeDefined()
    expect(table?.querySelector(`thead`)).toBeDefined()
    expect(table?.querySelector(`tbody`)).toBeDefined()

    // Check essential columns are present (with sort indicators)
    const header_texts = [...document.querySelectorAll(`th`)].map((h) =>
      h.textContent?.trim(),
    )
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
    expect(metric_order).toStrictEqual([...metric_order].sort((n1, n2) => n1 - n2))

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
      props: { col_filter: () => true, show_non_compliant: true },
    })
    await tick()

    const org_cell = doc_query(`td[data-col="Org"]`)
    const org_preview = doc_query(`td[data-col="Org"] .org-preview`)
    const headers = [...document.querySelectorAll(`th`)]
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
      props: { col_filter: () => true, show_non_compliant: true },
    })
    await tick()

    const cps_header = [...document.querySelectorAll(`th`)].find((header) =>
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
    const metadata_keys = new Set([`Training Set`, `Targets`, `date_added`, `Links`])
    // Labels displayed in table headers
    const metadata_labels = [`Training Set`, `Targets`, `Date Added`, `Links`]
    const col_filter = (_col: Label) => true // show all columns initially
    mount(MetricsTable, { target: document.body, props: { col_filter } })
    // Check metadata columns are visible initially
    let header_texts = [...document.querySelectorAll(`th`)].map((h) =>
      h.textContent?.replace(/\s*[↑↓]\s*$/, ``).trim(),
    )
    expect(header_texts).toStrictEqual(expect.arrayContaining(metadata_labels))

    // Create a new instance that hides metadata columns
    document.body.innerHTML = ``
    mount(MetricsTable, {
      target: document.body,
      props: {
        col_filter: (col: Label) => !metadata_keys.has(col.key ?? col.label),
        show_non_compliant: true,
      },
    })

    // Check metadata columns are hidden
    header_texts = [...document.querySelectorAll(`th`)].map((h) =>
      h.textContent?.replace(/\s*[↑↓]\s*$/, ``).trim(),
    )

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
      props: { col_filter, show_non_compliant: true },
    })

    const headers = document.querySelectorAll(`th`)
    const header_texts = [...headers].map((h) => h.textContent?.split(` `)[0])

    // Check hidden columns
    expect(header_texts).not.toContain(`F1`)
    expect(header_texts).not.toContain(`DAF`)

    // Check other columns still visible
    expect(header_texts).toContain(`Model`)
    expect(header_texts).toContain(`CPS`)
    expect(header_texts).toContain(`Prec`)
    expect(header_texts).toContain(`Acc`)
  })

  it(`filters energy-only models`, { timeout: 30_000 }, async () => {
    // First test with energy-only models hidden
    mount(MetricsTable, {
      target: document.body,
      props: { show_energy_only: false, show_non_compliant: true },
    })
    await tick()
    const rows_without_energy = document.querySelectorAll(`tbody tr`).length

    // Then test with energy-only models shown
    document.body.innerHTML = ``
    mount(MetricsTable, {
      target: document.body,
      props: { show_energy_only: true, show_non_compliant: true },
    })
    await tick()
    const rows_with_energy = document.querySelectorAll(`tbody tr`).length

    expect(rows_without_energy).toBe(visible_row_count({ energy: false }))
    expect(rows_with_energy).toBe(visible_row_count({ energy: true }))
  })

  it(`filters models based on model_filter prop`, () => {
    // First test: show no models
    const no_model_filter = (_model: ModelData) => false
    mount(MetricsTable, {
      target: document.body,
      props: {
        model_filter: no_model_filter,
        show_non_compliant: true,
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
        show_non_compliant: true,
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
        show_non_compliant: true,
      },
    })

    const filtered_rows = document.querySelectorAll(`tbody tr`)
    expect(filtered_rows).toHaveLength(
      visible_row_count({}, (model) => model.model_name.includes(`CHG`)),
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

    const headers = [...document.querySelectorAll(`th`)]
    expect(headers).toHaveLength(expected_headers.length)
    expect(headers.map((h) => h.textContent?.split(` `)[0])).toStrictEqual(
      expected_headers,
    )
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

      const headers = [...document.querySelectorAll(`th`)]
      expect(headers.map((h) => h.textContent?.split(` `)[0])).toStrictEqual(
        expected_headers,
      )

      const rows = document.querySelectorAll(`tbody tr`)
      rows.forEach((row) => {
        const model_cell = row.querySelector(`td`)
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
        show_non_compliant: true,
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

  it(`updates table when col_filter changes`, () => {
    // Test with only Model and F1 columns
    mount(MetricsTable, {
      target: document.body,
      props: {
        col_filter: (col: Label) => [`Model`, `F1`].includes(col.key ?? col.label),
        show_non_compliant: true,
      },
    })

    let headers = [...document.querySelectorAll(`th`)]
    expect(headers).toHaveLength(2)
    expect(headers.map((h) => h.textContent?.split(` `)[0])).toStrictEqual([
      `Model`,
      `F1`,
    ])

    // Create a new instance with Model, F1, and DAF columns
    document.body.innerHTML = ``
    mount(MetricsTable, {
      target: document.body,
      props: {
        col_filter: (col: Label) => [`Model`, `F1`, `DAF`].includes(col.key ?? col.label),
        show_non_compliant: true,
      },
    })

    headers = [...document.querySelectorAll(`th`)]
    expect(headers).toHaveLength(3)
    expect(headers.map((h) => h.textContent?.split(` `)[0])).toStrictEqual([
      `Model`,
      `F1`,
      `DAF`,
    ])
  })

  it(`renders table rows with default config`, async () => {
    mount(MetricsTable, {
      target: document.body,
      props: { show_non_compliant: true },
    })
    await tick()
    const rows = document.querySelectorAll(`tbody tr`)
    expect(rows).toHaveLength(visible_row_count())
  })

  describe(`Column Sorting`, () => {
    it(`sorts by Date Added chronologically, not alphabetically`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: {
          show_non_compliant: true,
          col_filter: (col: Label) => [`Model`, `Date Added`].includes(col.label),
        },
      })

      // Find Date Added column header
      const headers = [...document.querySelectorAll(`th`)]
      const date_header = headers.find((h) => h.textContent?.includes(`Date Added`))

      if (!date_header) {
        throw new Error(`Date Added column not found`)
      }

      // Click to sort by date
      date_header.click()
      await tick()

      // Get all date cells
      const date_cells = [
        ...document.querySelectorAll(`td[data-col="Date Added"] [data-sort-value]`),
      ]

      if (date_cells.length < 2) {
        throw new Error(`Not enough data for testing date sorting`)
      }

      // Get dates as timestamps from data-sort-value
      const timestamps = date_cells
        .map(parse_integer_sort_value)
        .filter((timestamp) => timestamp !== null)

      // rows must be fully sorted by timestamp (either direction) after the click
      const ascending = [...timestamps].toSorted((ts_1, ts_2) => ts_1 - ts_2)
      expect([ascending, ascending.toReversed()]).toContainEqual(timestamps)

      // Click again to toggle sort direction
      date_header.click()
      await tick()

      const reversed_timestamps = [
        ...document.querySelectorAll(`td[data-col="Date Added"] [data-sort-value]`),
      ]
        .map(parse_integer_sort_value)
        .filter((timestamp) => timestamp !== null)

      expect(reversed_timestamps).toStrictEqual(timestamps.toReversed())
    })

    it(`sorts numerically by training set size using data-sort-value`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: {
          show_non_compliant: true,
          col_filter: (col: Label) => [`Model`, `Training Set`].includes(col.label),
        },
      })

      // Find Training Set column header
      const headers = [...document.querySelectorAll(`th`)]
      const training_set_header = headers.find((h) =>
        h.textContent?.includes(`Training Set`),
      )

      if (!training_set_header) throw new Error(`Training Set column not found`)

      // Click to sort by training set size
      training_set_header.click()
      await tick()

      // Get training set sizes from data-sort-value
      const sizes = [
        ...document.querySelectorAll(`td[data-col="Training Set"] [data-sort-value]`),
      ]
        .map(parse_integer_sort_value)
        .filter((size) => size !== null)

      if (sizes.length < 2) {
        throw new Error(`Not enough data for testing training set sorting`)
      }

      // rows must be fully sorted by size (either direction) after the click
      const ascending = [...sizes].toSorted((size_1, size_2) => size_1 - size_2)
      expect([ascending, ascending.toReversed()]).toContainEqual(sizes)

      // Click again to toggle sort direction
      training_set_header.click()
      await tick()

      const reversed_sizes = [
        ...document.querySelectorAll(`td[data-col="Training Set"] [data-sort-value]`),
      ]
        .map(parse_integer_sort_value)
        .filter((size) => size !== null)

      expect(reversed_sizes).toStrictEqual(sizes.toReversed())
    })

    it(`sorts numerically by parameter count using data-sort-value`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: {
          show_non_compliant: true,
          col_filter: (col: Label) =>
            [`Model`, HYPERPARAMS.model_params.key].includes(col.key ?? col.label),
        },
      })

      // Find Params column header
      const headers = [...document.querySelectorAll(`th`)]
      const params_header = headers.find((h) => h.textContent?.includes(`Params`))

      if (!params_header) {
        throw new Error(`Params column not found`)
      }

      // Click to sort by parameter count
      params_header.click()
      await tick()

      // Get parameter counts from data-sort-value using the correct column label
      const param_counts = [
        ...document.querySelectorAll(
          `td[data-col="${HYPERPARAMS.model_params.label}"] [data-sort-value]`,
        ),
      ]
        .map(parse_integer_sort_value)
        .filter((count) => count !== null)

      if (param_counts.length < 2) {
        throw new Error(`Not enough data for testing parameter count sorting`)
      }

      // rows must be fully sorted by param count (either direction) after the click
      const ascending = [...param_counts].toSorted((cnt_1, cnt_2) => cnt_1 - cnt_2)
      expect([ascending, ascending.toReversed()]).toContainEqual(param_counts)

      // Click again to toggle sort direction
      params_header.click()
      await tick()

      const reversed_counts = [
        ...document.querySelectorAll(
          `td[data-col="${HYPERPARAMS.model_params.label}"] [data-sort-value]`,
        ),
      ]
        .map(parse_integer_sort_value)
        .filter((count) => count !== null)

      expect(reversed_counts).toStrictEqual(param_counts.toReversed())
    })

    it(`properly handles HTML content in cells without using it for data-sort-value`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: {
          show_non_compliant: true,
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
        props: { show_non_compliant: true, show_energy_only: true },
      },
      {
        test_name: `with filtered columns`,
        props: {
          show_non_compliant: true,
          col_filter: (col: Label) => [`Model`, `CPS`, `F1`].includes(col.label),
        },
      },
      {
        test_name: `with non-compliant models hidden`,
        props: { show_non_compliant: false, show_energy_only: true },
      },
    ])(
      `alphabetically sorts by Model name on $test_name header click`,
      { timeout: 30_000 }, // happy-dom renders of the full-column table are slow in CI
      async ({ props }) => {
        mount(MetricsTable, { target: document.body, props })

        // Find Model column header
        const headers = [...document.querySelectorAll(`th`)]
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
          visible_row_count({
            energy: `show_energy_only` in props && props.show_energy_only,
            non_compliant: props.show_non_compliant,
          }),
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
          show_non_compliant: true,
          col_filter: (col: Label) =>
            [`Model`, `CPS`, `Links`].includes(col.key ?? col.label),
        },
      })

      // Find CPS and Links column headers
      const headers = [...document.querySelectorAll(`th`)]
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
        props: { col_filter, show_non_compliant: true },
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
          show_non_compliant: true,
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
          show_non_compliant: true,
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
          show_non_compliant: true,
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
          props: { show_non_compliant: true },
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
    // diatomics each have Time and × Fastest columns, disambiguated by tooltip)
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
      `Time`, // ALL_METRICS (MD)
      `× Fastest`, // ALL_METRICS (MD)
      `CDS`, // DIATOMICS_METRICS
      `Time`, // DIATOMICS_METRICS
      `× Fastest`, // DIATOMICS_METRICS
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

    const header_elements = document.querySelectorAll(`thead th`)
    const actual_core_columns = [...header_elements].map((th) =>
      // Get text content, remove sort indicator (↑/↓) and any trailing spaces
      (th.textContent ?? ``).replace(/\s*[↑↓]\s*$/, ``).trim(),
    )

    // The default visible columns should stay intentionally curated: new default
    // columns must be added to expected_core_columns explicitly. Sorted comparison
    // ignores order but checks exact multiset (incl. duplicate Time/× Fastest labels).
    expect(actual_core_columns.toSorted()).toEqual(expected_core_columns.toSorted())

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
          props: { col_filter: () => true, show_non_compliant: true },
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
          props: { col_filter: () => true, show_non_compliant: true },
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
          props: { col_filter: () => true, show_non_compliant: true },
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
          props: { col_filter: () => true, show_non_compliant: true },
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

  describe(`regression tests for default values`, () => {
    // happy-dom renders of the full-column table are slow in CI
    it(
      `verifies critical default prop values to catch regressions`,
      {
        timeout: 30_000,
      },
      () => {
        mount(MetricsTable, {
          target: document.body,
          props: { col_filter: () => true, show_non_compliant: true },
        })

        // Verify table renders with data (filters allow content)
        const rows = document.querySelectorAll(`tbody tr`)
        expect(
          rows.length,
          `show_non_compliant=true & col_filter=true should show rows`,
        ).toBeGreaterThan(0)

        // Verify heatmap is enabled by default
        const table_controls = document.querySelector(`table-controls`)
        const heatmap_checkbox =
          table_controls?.querySelector<HTMLInputElement>(`input[type="checkbox"]`)
        expect(
          table_controls === null || heatmap_checkbox?.checked === true,
          `show_heatmap should default to true`,
        ).toBe(true)
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
          show_non_compliant: true,
        },
      })
      await tick()

      // After mounting, column_order should be initialized with ALL columns
      // (not just visible ones - the visible filter is separate)
      expect(state.column_order.length).toBeGreaterThan(10)
      expect(state.column_order).toContain(`Model`)
      expect(state.column_order).toContain(`F1`)
      expect(state.column_order).toContain(`DAF`)

      const visible_headers = [...document.querySelectorAll(`th`)].map(
        (h) => h.textContent?.split(` `)[0],
      )
      expect(visible_headers).toStrictEqual([`Model`, `F1`, `DAF`])
    })

    it.each([
      { columns: [`Model`, `F1`, `DAF`], name: `basic columns` },
      { columns: [`Model`, `F1`, `DAF`, `CPS`], name: `with CPS` },
    ])(`maintains Model column first with $name`, async ({ columns }) => {
      mount(MetricsTable, {
        target: document.body,
        props: {
          col_filter: (col: Label) => columns.includes(col.key ?? col.label),
          show_non_compliant: true,
        },
      })
      await tick()

      const headers = [...document.querySelectorAll(`th`)]
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
          show_non_compliant: true,
        },
      })
      await tick()

      const f1_idx = state.column_order.indexOf(`F1`)
      const daf_idx = state.column_order.indexOf(`DAF`)
      expect(f1_idx).toBeGreaterThanOrEqual(0)
      expect(daf_idx).toBeGreaterThanOrEqual(0)

      const headers = [...document.querySelectorAll(`th`)].map(
        (h) => h.textContent?.split(` `)[0],
      )
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
          show_non_compliant: true,
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

    it(`sets columns as draggable with correct attributes`, () => {
      mount(MetricsTable, {
        target: document.body,
        props: {
          col_filter: (col: Label) =>
            [`Model`, `F1`, `DAF`].includes(col.key ?? col.label),
          show_non_compliant: true,
        },
      })

      const headers = [...document.querySelectorAll(`thead tr:last-child th`)]
      expect(headers.length).toBeGreaterThan(0)
      headers.forEach((header) => {
        expect(header.getAttribute(`draggable`)).toBe(`true`)
        expect(header.getAttribute(`aria-dropeffect`)).toBe(`move`)
      })
    })

    it(`initializes without drag state classes`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: {
          col_filter: (col: Label) =>
            [`Model`, `F1`, `DAF`].includes(col.key ?? col.label),
          show_non_compliant: true,
        },
      })
      await tick()

      const headers = [...document.querySelectorAll<HTMLElement>(`th`)]

      // Initially no drag classes should be present
      headers.forEach((header) => {
        expect(header.classList.contains(`dragging`)).toBe(false)
        expect(header.classList.contains(`drag-over`)).toBe(false)
      })
    })
  })
})
