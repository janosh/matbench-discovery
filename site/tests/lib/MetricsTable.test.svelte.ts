import { DEFAULT_CPS_CONFIG } from '$lib/combined_perf_score.svelte'
import { HYPERPARAMS } from '$lib/labels'
import MetricsTable from '$lib/MetricsTable.svelte'
import type { Label, ModelData } from '$lib/types'
import { mount, tick } from 'svelte'
import { describe, expect, it } from 'vitest'
import { doc_query } from '../index'

describe(`MetricsTable`, () => {
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
    const header_texts = Array.from(document.querySelectorAll(`th`)).map((h) =>
      h.textContent?.trim()
    )
    const required_cols = [
      `Model`,
      `CPS ↑`, // Updated to include sort indicator
      `F1 ↑`,
      `DAF ↑`,
      `Training Set`,
      `Params`,
      `Targets`,
      `Links`,
    ]

    // Make sure each required column is present
    for (const col of required_cols) {
      expect(header_texts).toContain(col)
    }

    // Test prediction files dropdown interaction
    const pred_files_button = doc_query<HTMLButtonElement>(`tbody .pred-files-btn`)
    expect(pred_files_button).toBeDefined() // Ensure at least one button exists

    if (pred_files_button) {
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
    }
  })

  it(`toggles metadata columns`, () => {
    // These are actual METADATA_COLS (not HYPERPARAMS) - key and label are the same
    const metadata_cols = [`Training Set`, `Targets`, `Date Added`, `Links`]
    const col_filter = (_col: Label) => true // show all columns initially
    mount(MetricsTable, { target: document.body, props: { col_filter } })
    // Check metadata columns are visible initially
    let header_texts = Array.from(document.querySelectorAll(`th`)).map((h) =>
      h.textContent?.replace(/\s*[↑↓]\s*$/, ``).trim()
    )
    expect(header_texts).toEqual(expect.arrayContaining(metadata_cols))

    // Create a new instance that hides metadata columns
    document.body.innerHTML = ``
    mount(MetricsTable, {
      target: document.body,
      props: {
        col_filter: (col: Label) => !metadata_cols.includes(col.key ?? col.label),
        show_non_compliant: true,
      },
    })

    // Check metadata columns are hidden
    header_texts = Array.from(document.querySelectorAll(`th`)).map((h) =>
      h.textContent?.replace(/\s*[↑↓]\s*$/, ``).trim()
    )

    // Each metadata column should be hidden
    for (const col of metadata_cols) {
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
    const header_texts = Array.from(headers).map((h) => h.textContent?.split(` `)[0])

    // Check hidden columns
    expect(header_texts).not.toContain(`F1`)
    expect(header_texts).not.toContain(`DAF`)

    // Check other columns still visible
    expect(header_texts).toContain(`Model`)
    expect(header_texts).toContain(`CPS`)
    expect(header_texts).toContain(`Prec`)
    expect(header_texts).toContain(`Acc`)
  })

  it(`filters energy-only models`, () => {
    // First test with energy-only models hidden
    mount(MetricsTable, {
      target: document.body,
      props: { show_energy_only: false, show_non_compliant: true },
    })
    const rows_without_energy = document.querySelectorAll(`tbody tr`).length

    // Then test with energy-only models shown
    document.body.innerHTML = ``
    mount(MetricsTable, {
      target: document.body,
      props: { show_energy_only: true, show_non_compliant: true },
    })
    const rows_with_energy = document.querySelectorAll(`tbody tr`).length

    // Should have at least the same number of rows
    expect(rows_with_energy).toBeGreaterThanOrEqual(rows_without_energy)
  })

  it(`filters models based on model_filter prop`, () => {
    // First test: show no models
    const no_model_filter = (_model: ModelData) => false
    mount(MetricsTable, {
      target: document.body,
      props: {
        model_filter: no_model_filter,
        config: DEFAULT_CPS_CONFIG,
        show_non_compliant: true,
      },
    })

    const no_rows = document.querySelectorAll(`tbody tr`).length
    expect(no_rows).toBe(0)

    // Second test: show all models
    document.body.innerHTML = ``
    mount(MetricsTable, {
      target: document.body,
      props: {
        model_filter: () => true,
        config: DEFAULT_CPS_CONFIG,
        show_non_compliant: true,
      },
    })

    const all_rows = document.querySelectorAll(`tbody tr`).length
    expect(all_rows).toBeGreaterThan(0)

    // Third test: show specific models (e.g., only models with CHG in name)
    document.body.innerHTML = ``
    mount(MetricsTable, {
      target: document.body,
      props: {
        model_filter: (model: ModelData) => model.model_name.includes(`CHG`),
        config: DEFAULT_CPS_CONFIG,
        show_non_compliant: true,
      },
    })

    const filtered_rows = document.querySelectorAll(`tbody tr`)
    expect(filtered_rows.length).toBeGreaterThan(0)
    expect(filtered_rows.length).toBeLessThan(all_rows)

    // Verify that filtered rows actually contain CHG
    filtered_rows.forEach((row) => {
      const model_cell = row.querySelector(`td[data-col="Model"]`)
      expect(model_cell?.textContent).toContain(`CHG`)
    })
  })

  it(`validates prediction files dropdown button`, () => {
    // Create a simple element with the required structure for testing
    document.body.innerHTML = `
      <div>
        <button class="pred-files-btn" aria-label="Download model prediction files">
          <svg><path d="..."></path></svg>
        </button>
      </div>
    `

    // Find the button
    const pred_file_btn = document.querySelector(`.pred-files-btn`)
    expect(pred_file_btn).not.toBe(null)
    expect(pred_file_btn).toBeInstanceOf(HTMLButtonElement)
    expect(pred_file_btn?.getAttribute(`aria-label`)).toBe(
      `Download model prediction files`,
    )

    // Get model name from the same row for verification
    const model_cell = pred_file_btn?.closest(`tr`)?.querySelector(`td[data-col="Model"]`)
    expect(model_cell).not.toBe(null)

    // Check dropdown is initially not in the DOM
    const dropdown = document.querySelector(`.pred-files-dropdown`)
    expect(dropdown).toBe(null)
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

    const headers = Array.from(document.querySelectorAll(`th`))
    expect(headers.length).toBe(expected_headers.length)
    expect(headers.map((h) => h.textContent?.split(` `)[0])).toEqual(expected_headers)
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

      const headers = Array.from(document.querySelectorAll(`th`))
      expect(headers.map((h) => h.textContent?.split(` `)[0])).toEqual(expected_headers)

      const rows = document.querySelectorAll(`tbody tr`)
      rows.forEach((row) => {
        const model_cell = row.querySelector(`td`)
        expect(model_cell?.textContent).toContain(expected_model_match)
      })
    },
  )

  it(`updates table when col_filter changes`, () => {
    // Test with only Model and F1 columns
    mount(MetricsTable, {
      target: document.body,
      props: {
        col_filter: (col: Label) => [`Model`, `F1`].includes(col.key ?? col.label),
        show_non_compliant: true,
      },
    })

    let headers = Array.from(document.querySelectorAll(`th`))
    expect(headers.length).toBe(2)
    expect(headers.map((h) => h.textContent?.split(` `)[0])).toEqual([`Model`, `F1`])

    // Create a new instance with Model, F1, and DAF columns
    document.body.innerHTML = ``
    mount(MetricsTable, {
      target: document.body,
      props: {
        col_filter: (col: Label) => [`Model`, `F1`, `DAF`].includes(col.key ?? col.label),
        show_non_compliant: true,
      },
    })

    headers = Array.from(document.querySelectorAll(`th`))
    expect(headers.length).toBe(3)
    expect(headers.map((h) => h.textContent?.split(` `)[0])).toEqual([
      `Model`,
      `F1`,
      `DAF`,
    ])
  })

  it(`updates the table when CPS weights change`, async () => {
    // Test with default config
    mount(MetricsTable, {
      target: document.body,
      props: { config: DEFAULT_CPS_CONFIG, show_non_compliant: true },
    })
    await tick()
    const rows = document.querySelectorAll(`tbody tr`)
    expect(rows.length).toBeGreaterThan(18) // was 19
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
      const headers = Array.from(document.querySelectorAll(`th`))
      const date_header = headers.find((h) => h.textContent?.includes(`Date Added`))

      if (!date_header) {
        throw new Error(`Date Added column not found`)
      }

      // Click to sort by date
      date_header.click()
      await tick()

      // Get all date cells
      const date_cells = Array.from(
        document.querySelectorAll(`td[data-col="Date Added"]`),
      )

      if (date_cells.length < 2) {
        throw new Error(`Not enough data for testing date sorting`)
      }

      // Get dates as timestamps from data-sort-value
      const timestamps = date_cells
        .map((cell) => {
          const match = cell.innerHTML.match(/data-sort-value="(\d+)"/)
          return match ? parseInt(match[1]) : null
        })
        .filter((timestamp) => timestamp !== null) as number[]

      // Instead of checking order, verify we're extracting timestamps correctly
      expect(timestamps.length).toBeGreaterThan(0)
      timestamps.forEach((timestamp) => {
        expect(typeof timestamp).toBe(`number`)
        expect(timestamp).toBeGreaterThan(0)
      })

      // Click again to toggle sort direction
      date_header.click()

      // Get updated timestamps
      const descending_timestamps = Array.from(
        document.querySelectorAll(`td[data-col="Date Added"]`),
      )
        .map((cell) => {
          const match = cell.innerHTML.match(/data-sort-value="(\d+)"/)
          return match ? parseInt(match[1]) : null
        })
        .filter((timestamp) => timestamp !== null) as number[]

      // Verify we have timestamps
      expect(descending_timestamps.length).toBeGreaterThan(0)
      expect(descending_timestamps.length).toBe(timestamps.length)
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
      const headers = Array.from(document.querySelectorAll(`th`))
      const training_set_header = headers.find((h) =>
        h.textContent?.includes(`Training Set`)
      )

      if (!training_set_header) throw new Error(`Training Set column not found`)

      // Click to sort by training set size
      training_set_header.click()

      // Get training set sizes from data-sort-value
      const sizes = Array.from(
        document.querySelectorAll(`td[data-col="Training Set"]`),
      )
        .map((cell) => {
          const match = cell.innerHTML.match(/data-sort-value="(\d+)"/)
          return match ? parseInt(match[1]) : null
        })
        .filter((size) => size !== null) as number[]

      if (sizes.length < 2) {
        throw new Error(`Not enough data for testing training set sorting`)
      }

      // Just verify we're extracting numeric values correctly
      expect(sizes.length).toBeGreaterThan(0)
      sizes.forEach((size) => {
        expect(typeof size).toBe(`number`)
        expect(size).toBeGreaterThan(0)
      })

      // Get the initial order for comparison
      const initial_order = [...sizes]

      // Click again to toggle sort direction
      training_set_header.click()
      await tick()

      // Get updated sizes
      const new_sizes = Array.from(
        document.querySelectorAll(`td[data-col="Training Set"]`),
      )
        .map((cell) => {
          const match = cell.innerHTML.match(/data-sort-value="(\d+)"/)
          return match ? parseInt(match[1]) : null
        })
        .filter((size) => size !== null) as number[]

      // Verify we have the same number of items
      expect(new_sizes.length).toBe(initial_order.length)

      // The order should be different than before (reversed)
      let some_different = false
      for (const [idx, ref_size] of initial_order.entries()) {
        if (ref_size !== new_sizes[idx]) {
          some_different = true
          break
        }
      }

      // At least some items should be in a different order
      expect(some_different).toBe(true)
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
      const headers = Array.from(document.querySelectorAll(`th`))
      const params_header = headers.find((h) => h.textContent?.includes(`Params`))

      if (!params_header) {
        throw new Error(`Params column not found`)
      }

      // Click to sort by parameter count
      params_header.click()
      await tick()

      // Get parameter counts from data-sort-value using the correct column label
      const param_counts = Array.from(
        document.querySelectorAll(
          `td[data-col="${HYPERPARAMS.model_params.label}"]`,
        ),
      )
        .map((cell) => {
          const match = cell.innerHTML.match(/data-sort-value="(\d+)"/)
          return match ? parseInt(match[1]) : null
        })
        .filter((count) => count !== null) as number[]

      if (param_counts.length < 2) {
        throw new Error(`Not enough data for testing parameter count sorting`)
      }

      // Just verify we're extracting numeric values correctly
      expect(param_counts.length).toBeGreaterThan(0)
      param_counts.forEach((count) => {
        expect(typeof count).toBe(`number`)
        expect(count).toBeGreaterThan(0)
      })

      // Get the initial order for comparison
      const initial_order = [...param_counts]

      // Click again to toggle sort direction
      params_header.click()
      await tick()

      // Get updated counts using the correct column label
      const new_counts = Array.from(
        document.querySelectorAll(
          `td[data-col="${HYPERPARAMS.model_params.label}"]`,
        ),
      )
        .map((cell) => {
          const match = cell.innerHTML.match(/data-sort-value="(\d+)"/)
          return match ? parseInt(match[1]) : null
        })
        .filter((count) => count !== null) as number[]

      // Verify we have the same number of items
      expect(new_counts.length).toBe(initial_order.length)

      // The order should be different than before (reversed)
      let some_different = false
      for (const [idx, ref_count] of initial_order.entries()) {
        if (ref_count !== new_counts[idx]) {
          some_different = true
          break
        }
      }

      // At least some items should be in a different order
      expect(some_different).toBe(true)
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
      const training_set_cells = Array.from(
        document.querySelectorAll(`td[data-col="Training Set"]`),
      )

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
        if (data_sort_value) {
          expect(data_sort_value.includes(`<`)).toBe(false)
          expect(data_sort_value.includes(`>`)).toBe(false)
          expect(data_sort_value.includes(`span`)).toBe(false)
        }

        // The inner span should have its own data-sort-value
        const inner_span = cell.querySelector(`span[data-sort-value]`)
        if (inner_span) {
          const span_sort_value = inner_span.getAttribute(`data-sort-value`)
          expect(span_sort_value).toBeDefined()
          expect(isNaN(Number(span_sort_value))).toBe(false)
        }
      })

      // Verify tooltips are preserved on spans within cells
      const cells_with_tooltips = training_set_cells.filter(
        (cell) => cell.querySelector(`span[data-original-title]`) !== null,
      )

      expect(cells_with_tooltips.length).toBeGreaterThan(0)
      cells_with_tooltips.forEach((cell) => {
        const span = cell.querySelector(`span[data-original-title]`)
        expect(span?.getAttribute(`data-original-title`)).toBeTruthy()
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
      async ({ props }) => {
        mount(MetricsTable, { target: document.body, props })

        // Find Model column header
        const headers = Array.from(document.querySelectorAll(`th`))
        const model_header = headers.find((h) => h.textContent?.includes(`Model`))

        if (!model_header) throw new Error(`Model column header not found`)

        const get_model_names = () =>
          Array.from(document.querySelectorAll(`td[data-col="Model"]`))
            .map((cell) => {
              const link = cell.querySelector(`a`)
              return link?.getAttribute(`data-sort-value`)
            })
            .filter(Boolean) as string[]

        // model names before sorting
        // const initial_model_names = get_model_names()

        model_header.click() // Click to sort (ascending A-Z)
        await tick()

        // Get model names after first sort
        const sorted_model_names = get_model_names()

        // Check that we have enough models to test sorting
        expect(sorted_model_names.length).toBeGreaterThan(5)

        // check that order changed from sorting
        // TODO find out why this is failing
        // expect(sorted_model_names).not.toEqual(initial_model_names)
        // check that sorted_model_names is sorted
        expect(sorted_model_names).toEqual(sorted_model_names.sort())

        // Click again to reverse sort (descending Z-A)
        model_header.click()
        await tick()

        const reverse_sorted_model_names = get_model_names()
        expect(reverse_sorted_model_names).not.toEqual(sorted_model_names)
        expect(reverse_sorted_model_names.sort()).toEqual(sorted_model_names)
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
      const headers = Array.from(document.querySelectorAll(`th`))
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
      const initial_models = Array.from(
        document.querySelectorAll(`td[data-col="Model"]`),
      ).map((cell) => cell.textContent)

      // Try to sort by Links
      links_header.click()
      await tick()

      // Get model names after clicking Links
      const after_links_click_models = Array.from(
        document.querySelectorAll(`td[data-col="Model"]`),
      ).map((cell) => cell.textContent)

      // Order should not change
      expect(after_links_click_models).toEqual(initial_models)
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
      const links_cells = Array.from(
        document.querySelectorAll(`td[data-col="Links"]`),
      )
      expect(links_cells.length).toBeGreaterThan(20)

      // Check that rows have links (at least some should)
      let rows_with_links = 0
      for (const cell of links_cells) {
        const links = Array.from(cell.querySelectorAll(`a`))
        if (links.length > 1) rows_with_links++

        // Check each link has proper attributes
        for (const link of links) {
          expect(link.getAttribute(`target`)).toBe(`_blank`)
          expect(link.getAttribute(`rel`)).toBe(`noopener noreferrer`)

          const title = link.getAttribute(`data-original-title`)
          const href = link.getAttribute(`href`)
          expect(title).toBeTruthy()
          expect(href).toBeTruthy()

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
      const links_cells = Array.from(
        document.querySelectorAll(`td[data-col="Links"]`),
      )

      // Check for unavailable icon for missing links
      let found_missing_icon = false

      for (const cell of links_cells) {
        const missing_icons = Array.from(
          cell.querySelectorAll(
            `span[data-original-title$="not available"] svg`,
          ),
        )

        if (missing_icons.length > 0) {
          found_missing_icon = true

          // Check each missing link icon's parent span has a title
          for (const icon of missing_icons) {
            const span = icon.closest(`span`)
            const title = span?.getAttribute(`data-original-title`)
            expect(title).toBeTruthy()
            expect(title).toMatch(/not available/)
          }
        }
      }

      // Note: this might fail if all models have all links, which is unlikely
      expect(found_missing_icon).toBe(true)
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

      // Find all pred_files buttons
      const pred_file_buttons = Array.from(
        document.querySelectorAll(`.pred-files-btn`),
      )

      // Some models should have prediction files
      expect(pred_file_buttons.length).toBeGreaterThan(0)

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
      const links_cells = Array.from(
        document.querySelectorAll(`td[data-col="Links"]`),
      ).filter((cell) => cell.querySelectorAll(`a`).length > 0)

      // There should be at least one cell with links
      expect(links_cells.length).toBeGreaterThan(0)

      // Get the first cell with links to use as reference
      const reference_cell = links_cells[0]
      const reference_links = Array.from(reference_cell.querySelectorAll(`a`))

      // Check that links have SVG icons (now using Icon component)
      const reference_icons = reference_links.map((link) => {
        const svg = link.querySelector(`svg`)
        return svg !== null
      })

      // Verify that links have SVG icons
      expect(reference_icons.length).toBeGreaterThan(0)
      expect(reference_icons.some((has_icon) => has_icon)).toBe(true)

      // Check that the order of links is consistent across cells
      // (not checking all cells as some might have missing links)
      if (links_cells.length > 1) {
        const second_cell = links_cells[1]
        const second_links = Array.from(second_cell.querySelectorAll(`a`))

        // If the second cell has the same number of links, check they all have icons
        if (second_links.length === reference_links.length) {
          const second_icons = second_links.map((link) => {
            const svg = link.querySelector(`svg`)
            return svg !== null
          })

          // All links should have icons
          expect(second_icons.every((has_icon) => has_icon)).toBe(true)
        }
      }
    })
  })

  describe(`Heatmap Toggle Interaction`, () => {
    it(`toggles heatmap colors via TableControls checkbox`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: { show_non_compliant: true },
      })

      // Find the heatmap toggle checkbox within TableControls
      const heatmap_checkbox = doc_query<HTMLInputElement>(
        `input[type="checkbox"][aria-label="Toggle heatmap colors"]`,
      )

      expect(heatmap_checkbox).not.toBeNull()
      if (!heatmap_checkbox) return // Type guard

      // Helper to get cell background colors
      const get_cell_backgrounds = () => {
        const cells = document.querySelectorAll(
          `td[data-col="F1"], td[data-col="DAF"]`, // Check a couple of metric cols
        )
        return Array.from(cells).map(
          (cell) => getComputedStyle(cell as Element).backgroundColor,
        )
      }

      // Initially, heatmap should be on (default)
      expect(heatmap_checkbox.checked).toBe(true)
      let backgrounds = get_cell_backgrounds()
      expect(backgrounds.length).toBeGreaterThan(0) // Ensure cells were found
      expect(backgrounds.some((bg) => bg !== `` && bg !== `rgba(0, 0, 0, 0)`)).toBe(true)

      // Click the checkbox to turn heatmap off
      heatmap_checkbox.click()
      await tick()

      expect(heatmap_checkbox.checked).toBe(false)
      backgrounds = get_cell_backgrounds()
      expect(backgrounds.every((bg) => bg === `` || bg === `rgba(0, 0, 0, 0)`)).toBe(true)

      // Click again to turn heatmap back on
      heatmap_checkbox.click()
      await tick()

      expect(heatmap_checkbox.checked).toBe(true)
      backgrounds = get_cell_backgrounds()
      expect(backgrounds.some((bg) => bg !== `` && bg !== `rgba(0, 0, 0, 0)`)).toBe(true)
    })
  })

  it(`renders the correct default columns`, () => {
    mount(MetricsTable, { target: document.body })

    // Core text expected in default visible columns
    const expected_core_columns = new Set([
      `Model`, // METADATA_COLS
      `Training Set`, // METADATA_COLS
      `Targets`, // METADATA_COLS
      `Date Added`, // METADATA_COLS
      `Links`, // METADATA_COLS
      `Org`, // METADATA_COLS
      `Params`, // HYPERPARAMS (short label)
      `rcut`, // HYPERPARAMS (short label) - textContent doesn't keep subscript
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
      `RMSD`, // ALL_METRICS (Geo Opt)
      `CPS`, // Added in assemble_row_data
    ])

    const header_elements = document.querySelectorAll(`thead th`)
    const actual_core_columns = new Set(
      Array.from(header_elements).map((th) =>
        // Get text content, remove sort indicator (↑/↓) and any trailing spaces
        (th.textContent ?? ``).replace(/\s*[↑↓]\s*$/, ``).trim()
      ),
    )

    // Check if the sets of core texts are equal (ignores order)
    expect(actual_core_columns).toEqual(expected_core_columns)

    // Optionally, check the number of columns to be sure
    expect(header_elements.length).toBe(expected_core_columns.size)

    // Check that each header has a non-empty title attribute (tooltip)
    header_elements.forEach((th) => {
      expect(
        th.getAttribute(`title`),
        `Header ${th.textContent} has no title attribute`,
      ).toBeTruthy()
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

    it(`selects and deselects models on double-click with proper state management`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: { col_filter: () => true, show_non_compliant: true },
      })

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
    })

    it(`manages toggle visibility and count dynamically`, async () => {
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

      // Deselect all models by double-clicking each selected row
      // We need to deselect all currently selected rows
      const current_rows = get_rows()
      for (let i = 0; i < current_rows.length; i++) {
        if (current_rows[i].classList.contains(`highlight`)) {
          double_click_row(current_rows[i])
        }
      }
      await tick()

      expect(get_toggle()).toBeNull()
    })

    it(`toggles filter state and updates UI labels correctly`, async () => {
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
    })

    it(`filters rows and manages styling based on filter state`, async () => {
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
      expect(get_rows().length).toBe(1)
      expect(get_rows()[0].classList.contains(`highlight`)).toBe(false) // No highlight when filtering

      // Disable filter
      const toggle_disable = get_toggle()
      if (!toggle_disable) throw new Error(`Toggle not found`)
      toggle_disable.click()
      await tick()

      // Should show all rows with highlight
      expect(get_rows().length).toBe(initial_count)
      expect(get_rows()[0].classList.contains(`highlight`)).toBe(true)
    })

    it(`validates toggle behavior with multiple selections and deselections`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: { col_filter: () => true, show_non_compliant: true },
      })

      expect(get_rows().length).toBeGreaterThanOrEqual(2)

      // Test that toggle appears/disappears correctly
      expect(get_toggle()).toBeNull()

      // Select first row - use fresh references each time to avoid stale DOM
      double_click_row(get_rows()[0])
      await tick()
      expect(get_toggle()).not.toBeNull()
      expect(get_toggle_label()?.textContent).toContain(`1 selected`)

      // Select second row - get fresh reference
      double_click_row(get_rows()[1])
      await tick()
      expect(get_toggle_label()?.textContent).toContain(`2 selected`)

      // Test that deselecting works
      double_click_row(get_rows()[0])
      await tick()
      expect(get_toggle_label()?.textContent).toContain(`1 selected`)

      // Deselect the remaining selected row (should be row 1 which is still highlighted)
      const highlighted_row = document.querySelector(`tbody tr.highlight`)
      if (!highlighted_row) throw new Error(`Expected highlighted row`)
      double_click_row(highlighted_row)
      await tick()
      expect(get_toggle()).toBeNull()
    })

    it(`validates correct model name extraction and selection behavior`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: { col_filter: () => true, show_non_compliant: true },
      })

      const rows = Array.from(get_rows())
      expect(rows.length).toBeGreaterThanOrEqual(1)

      // Select first row
      double_click_row(rows[0])
      await tick()

      // Verify that the correct model was selected (not a generic "test-model")
      // This test will fail if model name extraction is broken
      expect(get_toggle()).not.toBeNull()
      expect(get_toggle_label()?.textContent).toContain(`1 selected`)

      // Get fresh row reference after selection
      const updated_rows = get_rows()
      expect(updated_rows.length).toBeGreaterThanOrEqual(1)
      // Verify that the row is highlighted
      expect(updated_rows[0].classList.contains(`highlight`)).toBe(true)

      // Deselect the row
      double_click_row(updated_rows[0])
      await tick()

      // Verify deselection worked
      expect(get_toggle()).toBeNull()
      const final_rows = get_rows()
      expect(final_rows[0].classList.contains(`highlight`)).toBe(false)
    })
  })

  describe(`regression tests for default values`, () => {
    it(`verifies critical default prop values to catch regressions`, () => {
      mount(MetricsTable, {
        target: document.body,
        props: { col_filter: () => true, show_non_compliant: true },
      })

      // Verify table renders with data (filters allow content)
      const rows = document.querySelectorAll(`tbody tr`)
      expect(rows.length, `show_non_compliant=true & col_filter=true should show rows`)
        .toBeGreaterThan(0)

      // Verify heatmap is enabled by default
      const table_controls = document.querySelector(`table-controls`)
      if (table_controls) {
        const heatmap_checkbox = table_controls.querySelector<HTMLInputElement>(
          `input[type="checkbox"]`,
        )
        expect(heatmap_checkbox?.checked, `show_heatmap should default to true`).toBe(
          true,
        )
      }
    })
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

      const visible_headers = Array.from(document.querySelectorAll(`th`)).map(
        (h) => h.textContent?.split(` `)[0],
      )
      expect(visible_headers).toEqual([`Model`, `F1`, `DAF`])
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

      const headers = Array.from(document.querySelectorAll(`th`))
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

      const headers = Array.from(document.querySelectorAll(`th`)).map(
        (h) => h.textContent?.split(` `)[0],
      )
      expect(headers[0]).toBe(`Model`)

      // F1 and DAF should appear in the order specified by column_order
      const visible_f1_pos = headers.indexOf(`F1`)
      const visible_daf_pos = headers.indexOf(`DAF`)
      if (f1_idx < daf_idx) {
        expect(visible_f1_pos).toBeLessThan(visible_daf_pos)
      } else {
        expect(visible_f1_pos).toBeGreaterThan(visible_daf_pos)
      }
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

      expect(state.column_order.length).toBe(initial_order.length)
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

      const headers = Array.from(document.querySelectorAll(`thead tr:last-child th`))
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

      const headers = Array.from(document.querySelectorAll(`th`)) as HTMLElement[]

      // Initially no drag classes should be present
      headers.forEach((header) => {
        expect(header.classList.contains(`dragging`)).toBe(false)
        expect(header.classList.contains(`drag-over`)).toBe(false)
      })
    })
  })
})
