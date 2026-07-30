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

    const checkbox = document.querySelector<HTMLInputElement>(`input#log`)
    const label = document.querySelector(`label[for="log"]`)

    expect(checkbox).toBeDefined()
    expect(label?.textContent).toContain(`Log color scale`)
  })

  it(`checkbox reflects log prop and toggles on click`, async () => {
    // Test default (false)
    mount(PtableHeatmap, {
      target: document.body,
      props: { heatmap_values: sample_values },
    })
    const checkbox = document.querySelector<HTMLInputElement>(`input#log`)
    expect(checkbox?.checked).toBe(false)

    // Test toggling
    checkbox?.click()
    await tick()
    expect(checkbox?.checked).toBe(true)

    checkbox?.click()
    await tick()
    expect(checkbox?.checked).toBe(false)

    // Test log=true prop
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

  // the colorbar range is [0, Math.max(0, ...values)]; without the 0 seed an empty
  // value set gives -Infinity and the colorbar renders no ticks at all
  it.each([
    { heatmap_values: {} as ElemCounts, desc: `empty`, ticks: [`0`] },
    {
      heatmap_values: { Fe: 100 } as ElemCounts,
      desc: `single element`,
      ticks: [`0`, `20`, `40`, `60`, `80`, `100`],
    },
    { heatmap_values: { Fe: 0, O: 0 } as ElemCounts, desc: `zero values`, ticks: [`0`] },
  ])(
    `renders finite colorbar ticks for $desc heatmap_values`,
    ({ heatmap_values, ticks }) => {
      mount(PtableHeatmap, { target: document.body, props: { heatmap_values } })

      const tick_labels = [...document.querySelectorAll(`.colorbar .tick-label`)].map(
        (label) => label.textContent?.trim(),
      )
      expect(tick_labels).toStrictEqual(ticks)
    },
  )
})
