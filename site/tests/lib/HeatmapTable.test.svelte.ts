import { HeatmapTable } from '$lib'
import type { Label } from '$lib/types'
import { mount, tick } from 'svelte'
import { describe, expect, it } from 'vitest'

describe(`HeatmapTable`, () => {
  const sample_data = [
    { Model: `Model A`, Score: 0.95, Value: 100 },
    { Model: `Model B`, Score: 0.85, Value: 200 },
    { Model: `Model C`, Score: 0.75, Value: 300 },
  ]

  const sample_columns: Label[] = [
    { key: `model`, label: `Model`, sticky: true, description: `` },
    { key: `score`, label: `Score`, better: `higher`, format: `.2f`, description: `` },
    { key: `value`, label: `Value`, better: `lower`, description: `` },
  ]

  it(`renders table with correct structure and handles hidden columns`, () => {
    const columns = [
      ...sample_columns,
      { key: `hidden`, label: `Hidden`, visible: false, description: `` },
    ]
    mount(HeatmapTable, {
      target: document.body,
      props: { data: sample_data, columns },
    })

    const headers = document.body.querySelectorAll(`th`)
    expect(headers).toHaveLength(3)
    expect(
      Array.from(headers).map((h) => h.textContent?.replace(/\s+/g, ` `).trim()),
    ).toEqual([`Model`, `Score ↑`, `Value ↓`])

    expect(document.body.querySelectorAll(`tbody tr`)).toHaveLength(3)
    expect(document.body.querySelectorAll(`td[data-col="Hidden"]`)).toHaveLength(0)
  })

  it(`handles empty data and filters undefined rows`, async () => {
    let data_with_empty = $state([{ Model: undefined, Score: undefined }, ...sample_data])

    mount(HeatmapTable, {
      target: document.body,
      props: { data: data_with_empty, columns: sample_columns },
    })

    expect(document.body.querySelectorAll(`tbody tr`)).toHaveLength(3)

    data_with_empty = []
    await tick()
    expect(document.body.querySelectorAll(`tbody tr`)).toHaveLength(3)
  })

  describe(`Sorting and Data Updates`, () => {
    it(`sorts correctly and handles missing values`, async () => {
      const data = [
        { Model: `A`, Score: undefined, Value: 100 },
        { Model: `B`, Score: 0.85, Value: undefined },
        { Model: `C`, Score: 0.75, Value: 300 },
      ]

      mount(HeatmapTable, {
        target: document.body,
        props: { data, columns: sample_columns },
      })

      // Test initial sort
      const value_header = document.body.querySelectorAll(`th`)[2]
      value_header.click()
      await tick()

      const values = Array.from(
        document.body.querySelectorAll(`td[data-col="Value"]`),
      ).map((cell) => cell.textContent?.trim())
      expect(values).toEqual([`100`, `300`, `n/a`])

      // Test sort direction toggle
      value_header.click()
      await tick()
      const reversed = Array.from(
        document.body.querySelectorAll(`td[data-col="Value"]`),
      ).map((cell) => cell.textContent?.trim())
      expect(reversed).toEqual([`300`, `100`, `n/a`])
    })

    it(`maintains sort state on data updates`, async () => {
      let data = $state(sample_data)
      mount(HeatmapTable, {
        target: document.body,
        props: { data, columns: sample_columns },
      })

      const score_header = document.body.querySelectorAll(`th`)[1]
      score_header.click() // Sort by Score
      await tick()

      data = [{ Model: `D`, Score: 0.65, Value: 400 }, ...sample_data]
      await tick()

      const scores = Array.from(
        document.body.querySelectorAll(`td[data-col="Score"]`),
      ).map((cell) => cell.textContent?.trim())
      expect(scores).toEqual([`0.95`, `0.85`, `0.75`])
    })

    it(`sorts date columns correctly`, () => {
      const dates = [
        { Date: `<span data-sort-value="1620950400000">2021-05-14</span>` },
        { Date: `<span data-sort-value="1684966800000">2023-05-25</span>` },
        { Date: `<span data-sort-value="1715089200000">2024-05-07</span>` },
      ]

      const date_columns: Label[] = [{ key: `date`, label: `Date`, description: `` }]

      mount(HeatmapTable, {
        target: document.body,
        props: { data: dates, columns: date_columns },
      })

      // Initial data should already be in order
      const initial_dates = Array.from(document.body.querySelectorAll(`td`)).map((cell) =>
        cell.textContent?.trim()
      )

      expect(initial_dates).toEqual([`2021-05-14`, `2023-05-25`, `2024-05-07`])
    })

    it(`sorts using data-sort-value attributes for numeric values`, () => {
      const formatted_data = [
        { Number: `<span data-sort-value="50">50</span>` },
        { Number: `<span data-sort-value="1000">1,000</span>` },
        { Number: `<span data-sort-value="10000">10,000</span>` },
      ]

      const columns: Label[] = [{ key: `number`, label: `Number`, description: `` }]

      mount(HeatmapTable, {
        target: document.body,
        props: { data: formatted_data, columns },
      })

      // Initial data
      const initial_numbers = Array.from(document.body.querySelectorAll(`td`)).map(
        (cell) => cell.textContent?.trim(),
      )
      expect(initial_numbers).toEqual([`50`, `1,000`, `10,000`])
    })

    it(`respects unsortable columns`, async () => {
      // Setup columns with an unsortable column
      const columns: Label[] = [
        { key: `name`, label: `Name`, sortable: true, description: `` },
        { key: `value`, label: `Value`, sortable: true, description: `` },
        { key: `actions`, label: `Actions`, sortable: false, description: `` },
      ]

      // Setup data with three sample entries
      const data = [
        { Name: `Alice`, Value: `100`, Actions: `View` },
        { Name: `Bob`, Value: `200`, Actions: `Edit` },
        { Name: `Charlie`, Value: `300`, Actions: `Delete` },
      ]

      mount(HeatmapTable, {
        target: document.body,
        props: { data, columns },
      })

      const headers = Array.from(document.body.querySelectorAll(`th`))
      const actions_header = headers[2]

      // Check initial values
      const initial_values = Array.from(
        document.body.querySelectorAll(`td[data-col="Value"]`),
      ).map((cell) => cell.textContent?.trim())

      expect(initial_values).toEqual([`100`, `200`, `300`])

      // Click the unsortable column - it should have no effect
      actions_header.click()
      await tick()

      // Capture values after clicking unsortable column
      const unchanged_values = Array.from(
        document.body.querySelectorAll(`td[data-col="Value"]`),
      ).map((cell) => cell.textContent?.trim())

      // Values should be unchanged since Actions is unsortable
      expect(unchanged_values).toEqual(initial_values)

      // Now try to sort by Value column
      headers[1].click()
      await tick()

      // Values should be sorted by Value column
      const post_sort_values = Array.from(
        document.body.querySelectorAll(`td[data-col="Value"]`),
      ).map((cell) => cell.textContent?.trim())

      // Check if values are sorted - the actual order depends on implementation
      expect(post_sort_values).not.toEqual(initial_values)
    })

    it(`sorts correctly with default initial sort settings`, () => {
      mount(HeatmapTable, {
        target: document.body,
        props: {
          data: sample_data,
          columns: sample_columns,
          initial_sort_column: `Score`,
          initial_sort_direction: `desc`,
        },
      })

      // Initial data should be sorted by Score in descending order
      const scores = Array.from(
        document.body.querySelectorAll(`td[data-col="Score"]`),
      ).map((cell) => cell.textContent?.trim())

      expect(scores).toEqual([`0.95`, `0.85`, `0.75`])
    })
  })

  it(`handles formatting and styles`, () => {
    const columns: Label[] = [
      { key: `num`, label: `Num`, format: `.1%`, description: `` },
      {
        key: `val`,
        label: `Val`,
        better: `higher`,
        color_scale: `interpolateViridis`,
        description: ``,
      },
    ]
    const data = [
      { Num: 0.123, Val: 0 },
      { Num: 1.234, Val: 100 },
    ]

    mount(HeatmapTable, {
      target: document.body,
      props: { data, columns },
    })

    // Check number formatting
    const num_cell = document.querySelector(`td[data-col="Num"]`)
    if (!num_cell) throw `Num cell not found`
    expect(num_cell.textContent?.trim()).toBe(`12.3%`)

    // Check that val cells have background colors
    const val_cells = document.querySelectorAll(`td[data-col="Val"]`)
    const backgrounds = Array.from(val_cells).map(
      (cell) => getComputedStyle(cell as Element).backgroundColor,
    )

    // Verify at least one background color is set (not empty)
    expect(backgrounds.some((bg) => bg !== `` && bg !== `rgba(0, 0, 0, 0)`)).toBe(true)
  })

  it(`applies different scale types for color mapping`, () => {
    const c1: Label = {
      key: `linear`,
      label: `Linear`,
      better: `higher`,
      color_scale: `interpolateViridis`,
      scale_type: `linear`,
      description: ``,
    }
    const c2: Label = {
      key: `log`,
      label: `Log`,
      better: `higher`,
      color_scale: `interpolateViridis`,
      scale_type: `log`,
      description: ``,
    }
    const data = [10, 100, 1000].map((val) => ({ [c1.label]: val, [c2.label]: val }))

    mount(HeatmapTable, { target: document.body, props: { data, columns: [c1, c2] } })

    // Get cells for both columns
    const linear_cells = document.querySelectorAll(`td[data-col="Linear"]`)
    const log_cells = document.querySelectorAll(`td[data-col="Log"]`)

    // Get background colors
    const linear_backgrounds = Array.from(linear_cells).map(
      (cell) => getComputedStyle(cell).backgroundColor,
    )
    const log_backgrounds = Array.from(log_cells).map(
      (cell) => getComputedStyle(cell).backgroundColor,
    )

    // Both types should have colors set
    expect(linear_backgrounds.every((bg) => bg !== `` && bg !== `rgba(0, 0, 0, 0)`)).toBe(
      true,
    )
    expect(log_backgrounds.every((bg) => bg !== `` && bg !== `rgba(0, 0, 0, 0)`)).toBe(
      true,
    )

    // The color distribution should be different between linear and log scale
    // In linear scale, 10->100->1000 should have increasingly spaced colors
    // In log scale, the color difference between 10->100 should be similar to 100->1000
    // Difficult to test precisely without mocking d3 scales, but we can check
    // there are differences between the two scale types.
    expect(linear_backgrounds).not.toEqual(log_backgrounds)
  })

  it(`handles accessibility features`, () => {
    mount(HeatmapTable, {
      target: document.body,
      props: {
        data: sample_data,
        columns: [{ key: `col`, label: `Col`, description: `Description`, sticky: true }],
        sort_hint: `Click to sort`,
      },
    })

    const header = document.body.querySelector(`th`)
    expect(header?.getAttribute(`title`) || header?.getAttribute(`data-title`)).toBe(
      `Description`,
    )
    expect(header?.querySelector(`[title="Click to sort"]`)).toBeDefined()
    expect(header?.classList.contains(`sticky-col`)).toBe(true)
  })

  it(`handles undefined and null values`, () => {
    const data = [{ Model: `Empty Model`, Score: undefined, Value: undefined }]

    mount(HeatmapTable, {
      target: document.body,
      props: { data, columns: sample_columns },
    })

    const cells = document.body.querySelectorAll(`td`)
    expect(cells).toHaveLength(3)
    expect(cells[0].textContent?.trim()).toBe(`Empty Model`)
    expect(cells[1].textContent?.trim()).toBe(`n/a`)
    expect(cells[2].textContent?.trim()).toBe(`n/a`)
  })

  it(`prevents HTML strings from being used as data-sort-value attributes`, () => {
    const html_data = [
      {
        Name: `Test Model`,
        HTML: `<span data-sort-value="100" title="This is a tooltip">100 units</span>`,
        Complex:
          `<span data-sort-value="3373529" title="Complex tooltip with multiple lines&#013;• Line item 1&#013;• Line item 2">3.37M <small>(details)</small> (<a href="https://example.com">Link</a>)</span>`,
      },
    ]

    const html_columns: Label[] = [
      { key: `name`, label: `Name`, description: `` },
      { key: `html`, label: `HTML`, description: `` },
      { key: `complex`, label: `Complex`, description: `` },
    ]

    mount(HeatmapTable, {
      target: document.body,
      props: { data: html_data, columns: html_columns },
    })

    // Get the cells with HTML content
    const html_cell = document.body.querySelector(`td[data-col="HTML"]`)
    const complex_cell = document.body.querySelector(`td[data-col="Complex"]`)

    // Verify cells exist and contain the expected HTML
    expect(html_cell).not.toBeNull()
    expect(complex_cell).not.toBeNull()

    // HTML should be rendered correctly
    expect(html_cell?.innerHTML).toContain(`<span data-sort-value="100"`)
    expect(complex_cell?.innerHTML).toContain(`<span data-sort-value="3373529"`)

    // The data-sort-value attribute on the td should not contain HTML
    const html_cell_sort_value = html_cell?.getAttribute(`data-sort-value`)
    const complex_cell_sort_value = complex_cell?.getAttribute(`data-sort-value`)

    // Either undefined (meaning HTML was detected and no sort value was set)
    // or not containing HTML tags
    if (html_cell_sort_value !== null) {
      expect(html_cell_sort_value?.includes(`<`)).toBe(false)
      expect(html_cell_sort_value?.includes(`>`)).toBe(false)
    }

    if (complex_cell_sort_value !== null) {
      expect(complex_cell_sort_value?.includes(`<`)).toBe(false)
      expect(complex_cell_sort_value?.includes(`>`)).toBe(false)
    }

    // Check that tooltips are present and accessible
    const tooltip_span = html_cell?.querySelector(`span[title]`)
    expect(tooltip_span).not.toBeNull()
    expect(tooltip_span?.getAttribute(`title`)).toBe(`This is a tooltip`)

    const complex_tooltip_span = complex_cell?.querySelector(`span[title]`)
    expect(complex_tooltip_span).not.toBeNull()
    expect(complex_tooltip_span?.getAttribute(`title`)).toContain(`Complex tooltip`)
  })

  describe(`Heatmap Toggle Functionality`, () => {
    it(`does not apply heatmap colors when show_heatmap is false`, () => {
      const columns: Label[] = [
        {
          key: `val`,
          label: `Val`,
          better: `higher`,
          color_scale: `interpolateViridis`,
          description: ``,
        },
      ]
      const data = [{ Val: 0 }, { Val: 100 }]

      mount(HeatmapTable, {
        target: document.body,
        props: { data, columns, show_heatmap: false }, // Disable heatmap
      })

      const val_cells = document.querySelectorAll(`td[data-col="Val"]`)
      const backgrounds = Array.from(val_cells).map(
        (cell) => getComputedStyle(cell as Element).backgroundColor,
      )

      // No background color should be applied when show_heatmap is false
      expect(backgrounds.every((bg) => bg === `` || bg === `rgba(0, 0, 0, 0)`)).toBe(true)
    })

    it(`applies heatmap colors when show_heatmap is true (default)`, () => {
      const columns: Label[] = [
        {
          key: `val`,
          label: `Val`,
          better: `higher`,
          color_scale: `interpolateViridis`,
          description: ``,
        },
      ]
      const data = [{ Val: 0 }, { Val: 100 }]

      mount(HeatmapTable, {
        target: document.body,
        props: { data, columns }, // show_heatmap is true by default
      })

      const val_cells = document.querySelectorAll(`td[data-col="Val"]`)
      const backgrounds = Array.from(val_cells).map(
        (cell) => getComputedStyle(cell as Element).backgroundColor,
      )

      // At least one background color should be set when show_heatmap is true
      expect(backgrounds.some((bg) => bg !== `` && bg !== `rgba(0, 0, 0, 0)`)).toBe(true)
    })
  })

  describe(`Column grouping`, () => {
    it(`correctly renders grouped columns`, () => {
      const grouped_columns: Label[] = [
        { key: `name`, label: `Name`, sticky: true, description: `` },
        { key: `val1_v`, label: `Value 1`, group: `Values`, description: `` },
        { key: `val2_v`, label: `Value 2`, group: `Values`, description: `` },
        { key: `met1`, label: `Metric 1`, group: `Metrics`, description: `` },
        { key: `met2`, label: `Metric 2`, group: `Metrics`, description: `` },
        { key: `val1_sv`, label: `Value 1`, group: `Second Values`, description: `` },
        { key: `val2_sv`, label: `Value 2`, group: `Second Values`, description: `` },
      ]

      const grouped_data = [
        {
          Name: `Item A`,
          'Value 1 (Values)': 10,
          'Value 2 (Values)': 20,
          'Metric 1': 30,
          'Metric 2': 40,
          'Value 1 (Second Values)': 50,
          'Value 2 (Second Values)': 60,
        },
      ]

      mount(HeatmapTable, {
        target: document.body,
        props: { data: grouped_data, columns: grouped_columns },
      })

      // Should have two rows in the header
      const header_rows = document.body.querySelectorAll(`thead tr`)
      expect(header_rows).toHaveLength(2)

      // First row should contain the group headers
      const group_headers = header_rows[0].querySelectorAll(`th`)
      expect(group_headers).toHaveLength(4) // Name (empty), Values, Metrics, Second Values

      // Get the text content of the group headers (excluding the empty one)
      const group_texts = Array.from(group_headers)
        .filter((th) => th.textContent?.trim())
        .map((th) => th.textContent?.trim())

      // Should have all three groups rendered
      expect(group_texts).toEqual([`Values`, `Metrics`, `Second Values`])

      // Check the group headers have correct colspan
      const values_header = Array.from(group_headers).find((th) =>
        th.textContent?.includes(`Values`)
      )
      const metrics_header = Array.from(group_headers).find((th) =>
        th.textContent?.includes(`Metrics`)
      )
      const second_values_header = Array.from(group_headers).find((th) =>
        th.textContent?.includes(`Second Values`)
      )

      expect(values_header?.getAttribute(`colspan`)).toBe(`2`) // Values spans 2 columns
      expect(metrics_header?.getAttribute(`colspan`)).toBe(`2`) // Metrics spans 2 columns
      expect(second_values_header?.getAttribute(`colspan`)).toBe(`2`) // Second Values spans 2 columns

      // Check column headers in second row
      const col_headers = header_rows[1].querySelectorAll(`th`)
      expect(col_headers).toHaveLength(7)

      // Column headers should have duplicate label names (Value 1, Value 2) rendered for each group
      expect(
        Array.from(col_headers).map((h) =>
          h.textContent?.trim().replace(/\s+|[↑↓]/g, ``)
        ),
      ).toEqual([`Name`, `Value1`, `Value2`, `Metric1`, `Metric2`, `Value1`, `Value2`])
    })

    it(`correctly handles mixed grouped and ungrouped columns`, () => {
      const mixed_columns: Label[] = [
        { key: `name`, label: `Name`, description: `` },
        { key: `regular`, label: `Regular`, description: `` },
        { key: `group1`, label: `Group 1`, group: `Grouped`, description: `` },
        { key: `group2`, label: `Group 2`, group: `Grouped`, description: `` },
        { key: `another`, label: `Another`, description: `` },
      ]

      const mixed_data = [
        { Name: `Test`, Regular: 1, 'Group 1': 2, 'Group 2': 3, Another: 4 },
      ]

      mount(HeatmapTable, {
        target: document.body,
        props: { data: mixed_data, columns: mixed_columns },
      })

      // Should have two rows in the header
      const header_rows = document.body.querySelectorAll(`thead tr`)
      expect(header_rows).toHaveLength(2)

      // Check the group header row
      const group_cells = header_rows[0].querySelectorAll(`th`)

      // There should be 4 cells - two empty (for Name and Regular), one for Grouped, and one empty for Another
      expect(group_cells).toHaveLength(4)

      // The Grouped cell should have colspan=2
      const grouped_header = Array.from(group_cells).find((c) =>
        c.textContent?.includes(`Grouped`)
      )
      expect(grouped_header?.getAttribute(`colspan`)).toBe(`2`)
    })
  })

  describe(`Style and CSS properties`, () => {
    it(`applies custom column styles`, () => {
      const styled_columns: Label[] = [
        {
          key: `col1`,
          label: `Col1`,
          style: `color: red; font-weight: lighter;`,
          description: ``,
        },
        { key: `col2`, label: `Col2`, description: `` },
      ]

      mount(HeatmapTable, {
        target: document.body,
        props: { data: [{ Col1: `a`, Col2: `b` }], columns: styled_columns },
      })

      const header = document.body.querySelector(`th`)
      expect(header?.getAttribute(`style`)).toContain(`color: red`)
      expect(header?.getAttribute(`style`)).toContain(`font-weight: lighter`)
      // Check that style is also applied to cells
      const cell = document.body.querySelector(`td[data-col="Col1"]`)
      expect(cell?.getAttribute(`style`)).toContain(`font-weight: lighter;`)
    })

    it(`applies row styles from data`, () => {
      const data_with_styles = [{ col: `value`, style: `background-color: yellow;` }]

      mount(HeatmapTable, {
        target: document.body,
        props: {
          data: data_with_styles,
          columns: [{ key: `col`, label: `col`, description: `` }],
        },
      })

      const row = document.body.querySelector(`tbody tr`)
      expect(row?.getAttribute(`style`)).toContain(`background-color: yellow`)
    })

    it(`applies container style from props`, () => {
      mount(HeatmapTable, {
        target: document.body,
        props: {
          data: [{ col: `value` }],
          columns: [{ key: `col`, label: `col`, description: `` }],
          style: `max-height: 200px; border: 1px solid blue;`,
        },
      })

      const container = document.body.querySelector(`.table-container`)
      expect(container?.getAttribute(`style`)).toContain(`max-height: 200px`)
      expect(container?.getAttribute(`style`)).toContain(`border: 1px solid blue`)
    })
  })
})
