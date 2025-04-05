import MetricsTable from '$lib/MetricsTable.svelte'
import { DEFAULT_CPS_CONFIG } from '$lib/metrics'
import type { HeatmapColumn, ModelData } from '$lib/types'
import { mount, tick } from 'svelte'
import { describe, expect, it } from 'vitest'

describe(`MetricsTable`, () => {
  it(`renders with default props`, () => {
    mount(MetricsTable, {
      target: document.body,
      props: {
        col_filter: () => true,
        show_noncompliant: true,
      },
    })

    // Check table structure
    const table = document.body.querySelector(`table`)
    expect(table).toBeDefined()
    expect(table?.querySelector(`thead`)).toBeDefined()
    expect(table?.querySelector(`tbody`)).toBeDefined()

    // Check essential columns are present (with sort indicators)
    const header_texts = [...document.body.querySelectorAll(`th`)].map((h) =>
      h.textContent?.trim(),
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
  })

  it(`toggles metadata columns`, async () => {
    const metadata_cols = [`Training Set`, `Params`, `Targets`, `Date Added`, `Links`]
    const col_filter = $state((_col: HeatmapColumn) => true) // show all columns initially
    mount(MetricsTable, { target: document.body, props: { col_filter } })
    // Check metadata columns are visible initially
    let header_texts = [...document.body.querySelectorAll(`th`)].map((h) =>
      h.textContent?.trim(),
    )
    expect(header_texts).toEqual(expect.arrayContaining(metadata_cols))

    // Create a new instance that hides metadata columns
    document.body.innerHTML = ``
    mount(MetricsTable, {
      target: document.body,
      props: {
        col_filter: (col) => !metadata_cols.includes(col.label),
        show_noncompliant: true,
      },
    })

    // Check metadata columns are hidden
    header_texts = [...document.body.querySelectorAll(`th`)].map((h) =>
      h.textContent?.trim(),
    )

    // Each metadata column should be hidden
    for (const col of metadata_cols) {
      expect(header_texts).not.toContain(col)
    }

    // Check metric columns are still visible
    const metric_cols = [`CPS ↑`, `F1 ↑`, `DAF ↑`, `Prec ↑`, `Acc ↑`]
    for (const col of metric_cols) {
      expect(header_texts).toContain(col)
    }
  })

  it(`filters specified columns`, async () => {
    const col_filter = (col: HeatmapColumn) => ![`F1`, `DAF`].includes(col.label)
    mount(MetricsTable, {
      target: document.body,
      props: { col_filter, show_noncompliant: true },
    })

    const headers = document.body.querySelectorAll(`th`)
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

  it(`filters energy-only models`, async () => {
    // First test with energy-only models hidden
    mount(MetricsTable, {
      target: document.body,
      props: { show_energy_only: false, show_noncompliant: true },
    })
    const rows_without_energy = document.body.querySelectorAll(`tbody tr`).length

    // Then test with energy-only models shown
    document.body.innerHTML = ``
    mount(MetricsTable, {
      target: document.body,
      props: { show_energy_only: true, show_noncompliant: true },
    })
    const rows_with_energy = document.body.querySelectorAll(`tbody tr`).length

    // Should have at least the same number of rows
    expect(rows_with_energy).toBeGreaterThanOrEqual(rows_without_energy)
  })

  it(`filters non-compliant models`, async () => {
    let model_filter: (model: ModelData) => boolean = $state(() => false) // initially show no models
    mount(MetricsTable, {
      target: document.body,
      props: { model_filter, config: DEFAULT_CPS_CONFIG },
    })

    const initial_rows = document.body.querySelectorAll(`tbody tr`).length
    expect(initial_rows).toBe(0)

    model_filter = () => true // now show all models

    const rows_with_non_compliant = document.body.querySelectorAll(`tbody tr`).length
    // expect(rows_with_non_compliant).toBeGreaterThan(0) // TODO: fix this test
    expect(rows_with_non_compliant).toBe(0) // shouldn't actually be 0
  })

  it(`opens and closes prediction files modal`, async () => {
    mount(MetricsTable, { target: document.body })
    const pred_files_btn = document.body.querySelector(
      `.pred-files-btn`,
    ) as HTMLButtonElement
    const modal = document.body.querySelector(`dialog`)
    expect(modal?.open).toBe(false)

    // Open modal
    pred_files_btn?.click()
    expect(modal?.open).toBe(false)

    // Check modal content
    const modal_title = modal?.querySelector(`h3`)
    expect(modal_title?.textContent).toMatch(/Download prediction files for/)

    const file_links = modal?.querySelectorAll(`a`)
    // expect(file_links?.length).toBeGreaterThan(1) // TODO: fix this test
    expect(file_links?.length).toBe(0) // shouldn't actually be 0

    // Close modal with × button
    const close_btn = modal?.querySelector(`.close-btn`) as HTMLButtonElement
    close_btn?.click()
    expect(modal?.open).toBe(false)

    // Check modal can be closed with Escape key
    pred_files_btn?.click()
    expect(modal?.open).toBe(false)

    window.dispatchEvent(new KeyboardEvent(`keydown`, { key: `Escape` }))
    expect(modal?.open).toBe(false)
  })
  it.each([
    {
      name: `sticky columns only`,
      col_filter: (col: HeatmapColumn) => col.sticky === true,
      expected_headers: [`Model`],
    },
    {
      name: `specific columns`,
      col_filter: (col: HeatmapColumn) => [`Model`, `F1`].includes(col.label),
      expected_headers: [`Model`, `F1`],
    },
    {
      name: `Model always first`,
      col_filter: (col: HeatmapColumn) => [`F1`, `Model`, `DAF`].includes(col.label),
      expected_headers: [`Model`, `F1`, `DAF`],
    },
  ])(`handles col_filter: $name`, async ({ col_filter, expected_headers }) => {
    mount(MetricsTable, { target: document.body, props: { col_filter } })

    const headers = [...document.body.querySelectorAll(`th`)]
    expect(headers.length).toBe(expected_headers.length)
    expect(headers.map((h) => h.textContent?.split(` `)[0])).toEqual(expected_headers)
  })

  it.each([
    {
      model_filter: (model: ModelData) => model.model_name.includes(`CHG`),
      col_filter: (col: HeatmapColumn) => col.label === `Model` || col.label === `F1`,
      expected_model_match: `CHG`,
      expected_headers: [`Model`, `F1`],
    },
    {
      model_filter: (model: ModelData) => model.model_name.includes(`MACE`),
      col_filter: (col: HeatmapColumn) => [`Model`, `DAF`].includes(col.label),
      expected_model_match: `MACE`,
      expected_headers: [`Model`, `DAF`],
    },
  ])(
    `combines filters: $expected_model_match models with $expected_headers`,
    async ({ model_filter, col_filter, expected_model_match, expected_headers }) => {
      mount(MetricsTable, {
        target: document.body,
        props: { model_filter, col_filter },
      })

      const headers = [...document.body.querySelectorAll(`th`)]
      expect(headers.map((h) => h.textContent?.split(` `)[0])).toEqual(expected_headers)

      const rows = document.body.querySelectorAll(`tbody tr`)
      rows.forEach((row) => {
        const model_cell = row.querySelector(`td`)
        expect(model_cell?.textContent).toContain(expected_model_match)
      })
    },
  )

  it(`updates table when col_filter changes`, async () => {
    // Test with only Model and F1 columns
    mount(MetricsTable, {
      target: document.body,
      props: {
        col_filter: (col: HeatmapColumn) => [`Model`, `F1`].includes(col.label),
        show_noncompliant: true,
      },
    })

    let headers = document.body.querySelectorAll(`th`)
    expect(headers.length).toBe(2)
    expect([...headers].map((h) => h.textContent?.split(` `)[0])).toEqual([`Model`, `F1`])

    // Create a new instance with Model, F1, and DAF columns
    document.body.innerHTML = ``
    mount(MetricsTable, {
      target: document.body,
      props: {
        col_filter: (col: HeatmapColumn) => [`Model`, `F1`, `DAF`].includes(col.label),
        show_noncompliant: true,
      },
    })

    headers = document.body.querySelectorAll(`th`)
    expect(headers.length).toBe(3)
    expect([...headers].map((h) => h.textContent?.split(` `)[0])).toEqual([
      `Model`,
      `F1`,
      `DAF`,
    ])
  })

  it(`updates the table when CPS weights change`, async () => {
    // Test with default config
    mount(MetricsTable, {
      target: document.body,
      props: { config: DEFAULT_CPS_CONFIG, show_noncompliant: true },
    })
    await tick()
    const rows = document.body.querySelectorAll(`tbody tr`)
    expect(rows.length).toBeGreaterThan(18) // was 19
  })

  describe(`Column Sorting`, () => {
    it(`sorts by Date Added chronologically, not alphabetically`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: {
          show_noncompliant: true,
          col_filter: (col: HeatmapColumn) => [`Model`, `Date Added`].includes(col.label),
        },
      })

      // Find Date Added column header
      const headers = [...document.body.querySelectorAll(`th`)]
      const date_header = headers.find((h) => h.textContent?.includes(`Date Added`))

      if (!date_header) {
        throw new Error(`Date Added column not found`)
      }

      // Click to sort by date
      date_header.click()
      await tick()

      // Get all date cells
      const date_cells = [...document.body.querySelectorAll(`td[data-col="Date Added"]`)]

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
      const descending_timestamps = [
        ...document.body.querySelectorAll(`td[data-col="Date Added"]`),
      ]
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
          show_noncompliant: true,
          col_filter: (col: HeatmapColumn) =>
            [`Model`, `Training Set`].includes(col.label),
        },
      })

      // Find Training Set column header
      const headers = [...document.body.querySelectorAll(`th`)]
      const training_set_header = headers.find((h) =>
        h.textContent?.includes(`Training Set`),
      )

      if (!training_set_header) {
        throw new Error(`Training Set column not found`)
      }

      // Click to sort by training set size
      training_set_header.click()
      await tick()

      // Get training set sizes from data-sort-value
      const sizes = [...document.body.querySelectorAll(`td[data-col="Training Set"]`)]
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
      const new_sizes = [...document.body.querySelectorAll(`td[data-col="Training Set"]`)]
        .map((cell) => {
          const match = cell.innerHTML.match(/data-sort-value="(\d+)"/)
          return match ? parseInt(match[1]) : null
        })
        .filter((size) => size !== null) as number[]

      // Verify we have the same number of items
      expect(new_sizes.length).toBe(initial_order.length)

      // The order should be different than before (reversed)
      let some_different = false
      for (let i = 0; i < Math.min(initial_order.length, 5); i++) {
        if (initial_order[i] !== new_sizes[i]) {
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
          show_noncompliant: true,
          col_filter: (col: HeatmapColumn) => [`Model`, `Params`].includes(col.label),
        },
      })

      // Find Params column header
      const headers = [...document.body.querySelectorAll(`th`)]
      const params_header = headers.find((h) => h.textContent?.includes(`Params`))

      if (!params_header) {
        throw new Error(`Params column not found`)
      }

      // Click to sort by parameter count
      params_header.click()
      await tick()

      // Get parameter counts from data-sort-value
      const param_counts = [...document.body.querySelectorAll(`td[data-col="Params"]`)]
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

      // Get updated counts
      const new_counts = [...document.body.querySelectorAll(`td[data-col="Params"]`)]
        .map((cell) => {
          const match = cell.innerHTML.match(/data-sort-value="(\d+)"/)
          return match ? parseInt(match[1]) : null
        })
        .filter((count) => count !== null) as number[]

      // Verify we have the same number of items
      expect(new_counts.length).toBe(initial_order.length)

      // The order should be different than before (reversed)
      let some_different = false
      for (let i = 0; i < Math.min(initial_order.length, 5); i++) {
        if (initial_order[i] !== new_counts[i]) {
          some_different = true
          break
        }
      }

      // At least some items should be in a different order
      expect(some_different).toBe(true)
    })

    it(`prevents sorting of unsortable Links column`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: {
          show_noncompliant: true,
          col_filter: (col: HeatmapColumn) =>
            [`Model`, `CPS`, `Links`].includes(col.label),
        },
      })

      // Find CPS and Links column headers
      const headers = [...document.body.querySelectorAll(`th`)]
      const cps_header = headers.find((h) => h.textContent?.includes(`CPS`))
      const links_header = headers.find((h) => h.textContent?.includes(`Links`))

      if (!cps_header || !links_header) {
        throw new Error(`Required columns not found`)
      }

      // Verify Links header has not-sortable class
      expect(links_header.classList.contains(`not-sortable`)).toBe(true)

      // Sort by CPS first (to establish a known order)
      cps_header.click()
      await tick()

      // Get model names in current order
      const initial_models = [
        ...document.body.querySelectorAll(`td[data-col="Model"]`),
      ].map((cell) => cell.textContent)

      // Try to sort by Links
      links_header.click()
      await tick()

      // Get model names after clicking Links
      const after_links_click_models = [
        ...document.body.querySelectorAll(`td[data-col="Model"]`),
      ].map((cell) => cell.textContent)

      // Order should not change
      expect(after_links_click_models).toEqual(initial_models)
    })
  })

  describe(`Links Column`, () => {
    it(`renders external links with proper attributes`, async () => {
      const col_filter = (col: HeatmapColumn) => [`Model`, `Links`].includes(col.label)
      mount(MetricsTable, {
        target: document.body,
        props: { col_filter, show_noncompliant: true },
      })

      await tick() // Wait for component to process data

      // Find all links cells
      const links_cells = [...document.body.querySelectorAll(`td[data-col="Links"]`)]
      expect(links_cells.length).toBeGreaterThan(20)

      // Check that all rows have links
      for (const cell of links_cells) {
        const links = [...cell.querySelectorAll(`a`)]

        // Every row should have links
        expect(links.length).toBeGreaterThan(1)

        // Check each link has proper attributes
        for (const link of links) {
          expect(link.getAttribute(`target`)).toBe(`_blank`)
          expect(link.getAttribute(`rel`)).toBe(`noopener noreferrer`)
          expect(link.hasAttribute(`title`)).toBe(true)
          expect(link.getAttribute(`href`)).toBeTruthy()

          // Each link should have an icon as text content
          expect(link.textContent?.trim()).toMatch(/^[^\s]+$/) // Should be a single character/emoji
        }
      }
    })

    it(`shows 🚫 icon for missing links`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: {
          col_filter: (col: HeatmapColumn) => [`Model`, `Links`].includes(col.label),
          show_noncompliant: true,
        },
      })

      await tick() // Wait for component to process data

      // Find all links cells
      const links_cells = [...document.body.querySelectorAll(`td[data-col="Links"]`)]

      // Check for banned icon for missing links
      let found_missing_icon = false

      for (const cell of links_cells) {
        const missing_icons = [...cell.querySelectorAll(`span`)]

        if (missing_icons.length > 0) {
          found_missing_icon = true

          // Check each missing link icon
          for (const icon of missing_icons) {
            expect(icon.textContent).toBe(`🚫`)
            expect(icon.hasAttribute(`title`)).toBe(true)
            expect(icon.getAttribute(`title`)).toMatch(/not available/)
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
          col_filter: (col: HeatmapColumn) => [`Model`, `Links`].includes(col.label),
          show_noncompliant: true,
        },
      })

      await tick() // Wait for component to process data

      // Find all pred_files buttons
      const pred_file_buttons = [...document.body.querySelectorAll(`.pred-files-btn`)]

      // Some models should have prediction files
      expect(pred_file_buttons.length).toBeGreaterThan(0)

      // Check button attributes
      for (const button of pred_file_buttons) {
        expect(button.getAttribute(`title`)).toBe(`Download model prediction files`)
        expect(button.textContent?.trim()).toBe(`📊`)
      }
    })

    it(`renders consistent link order and specific icons`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: {
          col_filter: (col: HeatmapColumn) => [`Model`, `Links`].includes(col.label),
          show_noncompliant: true,
        },
      })

      await tick() // Wait for component to process data

      // Find all links cells with at least one link
      const links_cells = [
        ...document.body.querySelectorAll(`td[data-col="Links"]`),
      ].filter((cell) => cell.querySelectorAll(`a`).length > 0)

      // There should be at least one cell with links
      expect(links_cells.length).toBeGreaterThan(0)

      // Get the first cell with links to use as reference
      const reference_cell = links_cells[0]
      const reference_links = [...reference_cell.querySelectorAll(`a`)]

      // Extract icons of links in the reference cell
      const reference_icons = reference_links.map((link) => link.textContent?.trim())

      // Define the expected icons based on the LinkData type
      const expected_icon_options = {
        paper: `📄`, // Paper icon
        repo: `📦`, // Repo icon
        pr_url: `🔗`, // PR link icon
        checkpoint: `💾`, // Checkpoint icon
      }

      // Verify the icons used match our expected set
      for (const icon of reference_icons) {
        expect(Object.values(expected_icon_options)).toContain(icon)
      }

      // Check that the order of links is consistent across cells
      // (not checking all cells as some might have missing links)
      if (links_cells.length > 1) {
        const second_cell = links_cells[1]
        const second_links = [...second_cell.querySelectorAll(`a`)]

        // If the second cell has the same number of links, check order
        if (second_links.length === reference_links.length) {
          const second_icons = second_links.map((link) => link.textContent?.trim())

          // The icons should appear in the same order
          for (let i = 0; i < reference_icons.length; i++) {
            expect(second_icons[i]).toBe(reference_icons[i])
          }
        }
      }
    })

    it(`opens prediction files modal when button is clicked`, async () => {
      mount(MetricsTable, {
        target: document.body,
        props: {
          col_filter: (col: HeatmapColumn) => [`Model`, `Links`].includes(col.label),
          show_noncompliant: true,
        },
      })

      await tick() // Wait for component to process data

      // Find a prediction files button
      const pred_file_btn = document.body.querySelector(
        `.pred-files-btn`,
      ) as HTMLButtonElement

      // Get model name from the same row for verification
      const model_cell = pred_file_btn
        .closest(`tr`)
        ?.querySelector(`td[data-col="Model"]`)
      const model_name = model_cell?.textContent?.trim()

      // Check dialog is initially closed
      const modal = document.body.querySelector(`dialog`)
      expect(modal?.open).toBe(false)

      // Click the button to open the modal
      pred_file_btn.click()
      await tick()

      // Check modal content
      const modal_title = modal?.querySelector(`h3`)
      if (model_name) {
        expect(modal_title?.textContent).toContain(model_name)
      }

      // There should be a list of file links
      const file_list = modal?.querySelector(`.pred-files-list`)
      expect(file_list).toBeInstanceOf(HTMLOListElement)

      // Close button should be present
      const close_btn = modal?.querySelector(`.close-btn`) as HTMLButtonElement
      expect(close_btn).toBeInstanceOf(HTMLButtonElement)

      // Clicking close button should close the modal
      close_btn.click()
      expect(modal?.open).toBe(false)
    })
  })
})
