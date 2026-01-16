import PtableHeatmap from '$lib/PtableHeatmap.svelte'
import type { ElementSymbol } from 'matterviz'
import { mount, tick } from 'svelte'
import { describe, expect, it } from 'vitest'

type ElemCounts = Record<ElementSymbol, number>
const sample_values = { H: 100, C: 250, O: 300, Fe: 150 } as ElemCounts

describe(`PtableHeatmap.svelte`, () => {
  it(`mounts and renders without errors`, () => {
    expect(() => {
      mount(PtableHeatmap, {
        target: document.body,
        props: { heatmap_values: sample_values },
      })
    }).not.toThrow()
    expect(document.body.innerHTML.length).toBeGreaterThan(0)
  })

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

    expect(component.snapshot.capture()).toEqual({
      color_scale: `interpolateViridis`,
      log: false,
    })

    component.snapshot.restore({ color_scale: `interpolatePlasma`, log: true })
    await tick()

    expect(component.snapshot.capture()).toEqual({
      color_scale: `interpolatePlasma`,
      log: true,
    })
  })

  it.each([
    { heatmap_values: {} as ElemCounts, desc: `empty` },
    { heatmap_values: { Fe: 100 } as ElemCounts, desc: `single element` },
    { heatmap_values: { Fe: 0, O: 0 } as ElemCounts, desc: `zero values` },
  ])(`mounts without error for $desc heatmap_values`, ({ heatmap_values }) => {
    expect(() => {
      mount(PtableHeatmap, { target: document.body, props: { heatmap_values } })
    }).not.toThrow()
  })
})
