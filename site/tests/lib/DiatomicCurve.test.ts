import DiatomicCurve from '$lib/plot/DiatomicCurve.svelte'
import { mount } from 'svelte'
import { describe, expect, it, vi } from 'vitest'
import { doc_query, get_scatter_plot_props } from '../index'

const plot_mocks = vi.hoisted(() => ({
  ScatterPlot: vi.fn(),
}))

vi.mock(`$lib`, () => ({ MODELS: [] }))
vi.mock(`matterviz`, () => ({ ScatterPlot: plot_mocks.ScatterPlot }))

describe(`DiatomicCurve`, () => {
  it(`passes filtered and shifted curve data to ScatterPlot`, () => {
    mount(DiatomicCurve, {
      target: document.body,
      props: {
        formula: `H2`,
        class: `custom-class`,
        curves: [
          {
            model_key: `unknown-model`,
            distances: [0.1, 0.2, 1.0, 6.0, 6.1],
            energies: [100, 5, 8, 3, 100],
            color: `#123456`,
          },
        ],
      },
    })

    expect(doc_query(`h3`).childNodes[0]?.textContent?.trim()).toBe(`H2`)
    expect(doc_query(`.plot`).classList.contains(`custom-class`)).toBe(true)
    expect(plot_mocks.ScatterPlot).toHaveBeenCalledTimes(1)
    expect(get_scatter_plot_props(plot_mocks.ScatterPlot)).toMatchObject({
      series: [
        {
          x: [0.2, 1.0, 6.0],
          y: [2, 5, 0],
          markers: `line+points`,
          metadata: { model_key: `unknown-model`, model_label: `unknown-model` },
          point_style: { fill: `#123456`, radius: 1.5, stroke_width: 0 },
          line_style: { stroke: `#123456` },
        },
      ],
      x_axis: { label: `Distance (Å)`, format: `.1f`, range: [0.2, 6] },
      y_axis: { label: `Energy (eV)`, format: `.2f`, range: [-8, 20] },
      legend: null,
    })
  })
})
