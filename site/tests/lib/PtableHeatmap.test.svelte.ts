import PtableHeatmap from '$lib/PtableHeatmap.svelte'
import type { ElementSymbol } from 'matterviz'
import { tick } from 'svelte'
import { describe, expect, it } from 'vitest'
import { mount } from '../index'

type ElemCounts = Record<ElementSymbol, number>
const sample_values = { H: 100, C: 250, O: 300, Fe: 150 } as ElemCounts

describe(`PtableHeatmap.svelte`, () => {
  it(`renders log scale checkbox with label`, () => {
    mount(PtableHeatmap, {
      target: document.body,
      props: { heatmap_values: sample_values },
    })

    const checkbox = document.querySelector(`input#log`)
    expect(checkbox).toBeInstanceOf(HTMLInputElement)
    expect(document.querySelector(`label[for="log"]`)?.textContent).toContain(
      `Log color scale`,
    )
  })

  it(`checkbox reflects log prop and toggles on click`, async () => {
    mount(PtableHeatmap, {
      target: document.body,
      props: { heatmap_values: sample_values },
    })
    const checkbox = document.querySelector<HTMLInputElement>(`input#log`)
    expect(checkbox?.checked).toBe(false)

    checkbox?.click()
    await tick()
    expect(checkbox?.checked).toBe(true)

    checkbox?.click()
    await tick()
    expect(checkbox?.checked).toBe(false)

    document.body.innerHTML = ``
    mount(PtableHeatmap, {
      target: document.body,
      props: { heatmap_values: sample_values, log: true },
    })
    expect(document.querySelector<HTMLInputElement>(`input#log`)?.checked).toBe(true)
  })

  it(`snapshot capture and restore work correctly`, async () => {
    const component = mount(PtableHeatmap, {
      target: document.body,
      props: {
        heatmap_values: sample_values,
        color_scale: `interpolateViridis`,
        log: false,
      },
    })

    expect(component.snapshot.capture()).toStrictEqual({
      color_scale: `interpolateViridis`,
      log: false,
    })

    component.snapshot.restore({ color_scale: `interpolatePlasma`, log: true })
    await tick()

    expect(component.snapshot.capture()).toStrictEqual({
      color_scale: `interpolatePlasma`,
      log: true,
    })
  })

  // colorbar range is [0, Math.max(0, ...values)]; without the 0 seed, empty → -Infinity
  it.each([
    { vals: {} as ElemCounts, max: 0 },
    { vals: { Fe: 100 } as ElemCounts, max: 100 },
    { vals: { Fe: 0, O: 0 } as ElemCounts, max: 0 },
  ])(`colorbar max tick is $max`, ({ vals, max }) => {
    mount(PtableHeatmap, { target: document.body, props: { heatmap_values: vals } })
    const ticks = [...document.querySelectorAll(`.colorbar .tick-label`)].map((label) =>
      Number(label.textContent),
    )
    expect(Math.max(...ticks)).toBe(max)
  })
})
