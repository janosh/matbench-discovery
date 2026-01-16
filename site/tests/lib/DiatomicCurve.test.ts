import { DiatomicCurve } from '$lib'
import { mount } from 'svelte'
import { describe, expect, it, vi } from 'vitest'
import { doc_query } from '../index'

vi.mock(`matterviz`, () => ({ ScatterPlot: vi.fn() }))

const create_mock_curve = (overrides: Partial<{
  model_key: string
  distances: number[]
  energies: number[]
  color: string
}> = {}) => ({
  model_key: `test-model`,
  distances: [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0],
  energies: [-2.0, -5.0, -6.0, -4.0, -2.0, -1.0, -0.2, -0.05, 0.0],
  color: `#4285F4`,
  ...overrides,
})

describe(`DiatomicCurve`, () => {
  it(`renders formula in heading`, () => {
    mount(DiatomicCurve, {
      target: document.body,
      props: { formula: `H2`, curves: [create_mock_curve()] },
    })
    expect(doc_query<HTMLHeadingElement>(`h3`).textContent).toBe(`H2`)
  })

  it(`applies custom class via rest props`, () => {
    mount(DiatomicCurve, {
      target: document.body,
      props: { formula: `H2`, curves: [create_mock_curve()], class: `custom-class` },
    })
    expect(doc_query<HTMLDivElement>(`.plot`).classList.contains(`custom-class`)).toBe(
      true,
    )
  })
})
