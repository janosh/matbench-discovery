import { HeatmapTable } from '$lib'
import type { HeatmapColumn } from '$lib/types'
import { mount, tick } from 'svelte'
import { describe, expect, it } from 'vitest'

describe(`HeatmapTable`, () => {
  const sample_data = [
    { Model: `Model A`, Score: 0.95, Value: 100 },
    { Model: `Model B`, Score: 0.85, Value: 200 },
    { Model: `Model C`, Score: 0.75, Value: 300 },
  ]

  const sample_columns: HeatmapColumn[] = [
    { label: `Model`, sticky: true },
    { label: `Score`, better: `higher`, format: `.2f` },
    { label: `Value`, better: `lower` },
  ]

  it(`renders table with correct structure and handles hidden columns`, () => {
    const columns = [...sample_columns, { label: `Hidden`, hidden: true }]
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
    const data_with_empty = $state([
      { Model: undefined, Score: undefined },
      ...sample_data,
    ])

    mount(HeatmapTable, {
      target: document.body,
      props: { data: data_with_empty, columns: sample_columns },
    })

    expect(document.body.querySelectorAll(`tbody tr`)).toHaveLength(3)

    data_with_empty = []
    await tick()
    expect(document.body.querySelectorAll(`tbody tr`)).toHaveLength(0)
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
      const component = mount(HeatmapTable, {
        target: document.body,
        props: { data: sample_data, columns: sample_columns },
      })

      const score_header = document.body.querySelectorAll(`th`)[1]
      score_header.click() // Sort by Score
      await tick()

      component.$set({ data: [{ Model: `D`, Score: 0.65 }, ...sample_data] })
      await tick()

      const scores = Array.from(
        document.body.querySelectorAll(`td[data-col="Score"]`),
      ).map((cell) => cell.textContent?.trim())
      expect(scores).toEqual([`0.65`, `0.95`, `0.85`, `0.75`])
    })
  })

  it(`handles formatting and styles`, () => {
    const columns: HeatmapColumn[] = [
      { label: `Num`, format: `.1%` },
      { label: `Val`, better: `higher` },
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
    expect(getComputedStyle(num_cell).backgroundColor).toBe(`rgb(68, 1, 84)`)

    // Check value cells have different background colors
    const val_cells = document.querySelectorAll(`td[data-col="Val"]`)
    const backgrounds = Array.from(val_cells).map(
      (cell) => getComputedStyle(cell).backgroundColor,
    )
    expect(backgrounds).toEqual([`rgb(68, 1, 84)`, `rgb(253, 231, 37)`])
  })

  it(`handles accessibility features`, () => {
    mount(HeatmapTable, {
      target: document.body,
      props: {
        data: sample_data,
        columns: [{ label: `Col`, tooltip: `Description`, sticky: true }],
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
})
