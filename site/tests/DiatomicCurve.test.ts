import { DiatomicCurve } from '$lib'
import { mount } from 'svelte'
import { beforeEach, describe, expect, test } from 'vitest'
import { doc_query } from '.'

describe(`DiatomicCurve`, () => {
  // Add container with dimensions to body before each test
  const container_style = `width: 800px; height: 600px;`
  beforeEach(() => {
    const container = document.createElement(`div`)
    container.setAttribute(`style`, container_style)
    document.body.appendChild(container)
  })

  const sample_curves = [
    {
      model_name: `mace-mp-0-medium`,
      distances: [0.5, 1.0, 1.5, 2.0, 2.5],
      energies: [-5, -7, -3, -1, 0],
      color: `#ff0000`,
    },
    {
      model_name: `mace-mpa-0-medium`,
      distances: [0.5, 1.0, 1.5, 2.0, 2.5],
      energies: [-4, -6, -2, -0.5, 0],
      color: `#00ff00`,
    },
  ]

  test(`renders with basic props`, async () => {
    mount(DiatomicCurve, {
      target: document.body,
      props: { formula: `H₂`, curves: sample_curves },
    })

    const title = doc_query(`h3`)
    expect(title.textContent).toBe(`H₂`)
  })

  test(`applies correct styling`, async () => {
    const style = `background: #f0f0f0;`
    mount(DiatomicCurve, {
      target: document.body,
      props: { formula: `H₂`, curves: sample_curves, style },
    })

    const plot = doc_query(`.plot`)
    expect(plot.getAttribute(`style`)).toContain(style)
  })
})
