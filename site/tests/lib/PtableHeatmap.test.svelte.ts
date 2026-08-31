import PtableHeatmap from '$lib/PtableHeatmap.svelte'
import { format_num, type ElementSymbol } from 'matterviz'
import { tick } from 'svelte'
import { describe, expect, it } from 'vitest'
import { checkbox_for, mount } from '../index'

type ElemCounts = Record<ElementSymbol, number>
const sample_values = { H: 100, C: 250, O: 300, Fe: 150 } as ElemCounts

const colorbar_ticks = (): number[] =>
  [...document.querySelectorAll(`.colorbar .tick-label`)].map((label) =>
    Number(label.textContent),
  )

describe(`PtableHeatmap.svelte`, () => {
  it(`renders labeled log checkbox that reflects log prop and toggles on click`, async () => {
    mount(PtableHeatmap, {
      target: document.body,
      props: { heatmap_values: sample_values },
    })
    const checkbox = checkbox_for(`Log color scale`)
    expect(checkbox.checked).toBe(false)

    checkbox.click()
    await tick()
    expect(checkbox.checked).toBe(true)

    checkbox.click()
    await tick()
    expect(checkbox.checked).toBe(false)

    document.body.innerHTML = ``
    mount(PtableHeatmap, {
      target: document.body,
      props: { heatmap_values: sample_values, log: true },
    })
    expect(checkbox_for(`Log color scale`).checked).toBe(true)
  })

  // /data mounts three heatmaps; a global id would be duplicated (and label[for] would
  // point every label at the first checkbox)
  it(`nests the log checkbox in its label without a global id`, () => {
    mount(PtableHeatmap, {
      target: document.body,
      props: { heatmap_values: sample_values },
    })
    expect(document.querySelector(`input#log`)).toBeNull()
    expect(document.querySelector(`label[for]`)).toBeNull()
    expect(checkbox_for(`Log color scale`).type).toBe(`checkbox`)
  })

  // colorbar range spans the finite data extent; without the 0 fallback, empty → -Infinity
  it.each([
    { vals: {} as ElemCounts, min: 0, max: 0 },
    { vals: { Fe: 100 } as ElemCounts, min: 100, max: 100 },
    { vals: { Fe: 0, O: 0 } as ElemCounts, min: 0, max: 0 },
    { vals: sample_values, min: 100, max: 300 },
  ])(`colorbar ticks span [$min, $max]`, ({ vals, min, max }) => {
    mount(PtableHeatmap, { target: document.body, props: { heatmap_values: vals } })
    const ticks = colorbar_ticks()
    expect(Math.min(...ticks)).toBe(min)
    expect(Math.max(...ticks)).toBe(max)
  })

  // in log mode the legend must be log-spaced (like the tiles) and skip non-positive values
  it(`log mode draws a log-spaced legend over the positive values`, async () => {
    const heatmap_values = { H: 0, C: 1, O: 10_000, Fe: -5 } as ElemCounts
    mount(PtableHeatmap, { target: document.body, props: { heatmap_values, log: true } })
    await tick()
    // 5 ticks evenly spaced in log space between 1 and 1e4 → decades (SI-formatted labels)
    const tick_texts = [...document.querySelectorAll(`.colorbar .tick-label`)].map(
      (label) => label.textContent,
    )
    expect(tick_texts).toEqual([1, 10, 100, 1000, 10_000].map((val) => format_num(val)))

    document.body.innerHTML = ``
    mount(PtableHeatmap, { target: document.body, props: { heatmap_values, log: false } })
    // linear ticks span the full finite extent, so the (niced) low end is negative
    const linear_ticks = [...document.querySelectorAll(`.colorbar .tick-label`)].map(
      (label) => label.textContent,
    )
    expect(linear_ticks[0]?.startsWith(`−`)).toBe(true)
    expect(linear_ticks).not.toEqual(tick_texts)
  })
})
