import { DiatomicCurve } from '$lib'
import { mount } from 'svelte'
import { beforeEach, describe, expect, it } from 'vitest'

describe(`DiatomicCurve`, () => {
  // Add container with dimensions to body before each test
  beforeEach(() => {
    document.body
      .appendChild(document.createElement(`div`))
      .setAttribute(`style`, `width: 800px; height: 600px;`)
  })

  const sample_curves = [
    {
      model_key: `mace-mp-0-medium`,
      distances: [0.5, 1.0, 1.5, 2.0, 2.5],
      energies: [-5, -7, -3, -1, 0],
      color: `#ff0000`,
    },
    {
      model_key: `mace-mpa-0-medium`,
      distances: [0.5, 1.0, 1.5, 2.0, 2.5],
      energies: [-4, -6, -2, -0.5, 0],
      color: `#00ff00`,
    },
  ]

  it(`renders with basic props`, () => {
    mount(DiatomicCurve, {
      target: document.body,
      props: { formula: `H₂`, curves: sample_curves },
    })

    // Verify component renders with expected content
    expect(document.querySelector(`h3, .title, .formula`)).toBeDefined()
    expect(document.body.textContent).toContain(`H₂`)
  })

  it(`applies correct styling`, () => {
    const style = `background: rgb(240, 240, 240); color: rgb(0, 0, 0);`
    mount(DiatomicCurve, {
      target: document.body,
      props: { formula: `H₂`, curves: sample_curves, style },
    })

    const curve_container = document.querySelector(`.plot`)
    expect(curve_container).toBeDefined()
    expect(curve_container?.getAttribute(`style`)).toContain(style)
  })
})
