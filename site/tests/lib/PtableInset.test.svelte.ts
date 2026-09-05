import PtableInset from '$lib/PtableInset.svelte'
import type { ChemicalElement } from 'matterviz'
import { describe, expect, it } from 'vitest'
import { mount } from '../index'

const mock_Fe = { symbol: `Fe`, name: `Iron`, number: 26 } as ChemicalElement
const mock_H = { symbol: `H`, name: `Hydrogen`, number: 1 } as ChemicalElement

describe(`PtableInset.svelte`, () => {
  it.each([
    { source: `record`, element: mock_Fe, elem_counts: { Fe: 150, O: 50 }, count: 150 },
    {
      source: `array indexed by atomic number`,
      element: mock_H,
      elem_counts: [200, ...Array<number>(119).fill(0)],
      count: 200,
    },
  ])(`renders element name and count from $source`, ({ element, elem_counts, count }) => {
    mount(PtableInset, { target: document.body, props: { element, elem_counts } })

    const strong = document.querySelector(`strong`)
    expect(strong?.textContent).toContain(element.name)
    expect(strong?.textContent).toContain(String(count))
  })

  it.each([undefined, false])(`shows percentage when show_percent=%s`, (show_percent) => {
    mount(PtableInset, {
      target: document.body,
      props: { element: mock_Fe, elem_counts: { Fe: 50, O: 50 }, show_percent },
    })
    expect(document.querySelector(`strong`)?.textContent?.includes(`%`)).toBe(
      show_percent !== false,
    )
  })

  it(`displays unit and renders HTML in unit prop`, () => {
    mount(PtableInset, {
      target: document.body,
      props: {
        element: mock_Fe,
        elem_counts: { Fe: 100 },
        unit: `<sub>2</sub>`,
        show_percent: false,
      },
    })

    const strong = document.querySelector(`strong`)
    expect(strong?.innerHTML).toContain(`<sub>2</sub>`)
  })

  it(`forwards class and style props`, () => {
    mount(PtableInset, {
      target: document.body,
      props: {
        element: mock_Fe,
        elem_counts: { Fe: 100 },
        class: `custom-class`,
        style: `color: red;`,
      },
    })

    const strong = document.querySelector(`strong`)
    expect(strong?.classList.contains(`custom-class`)).toBe(true)
    expect(strong?.getAttribute(`style`)).toContain(`color: red`)
  })

  it(`renders nothing when element has no name`, () => {
    mount(PtableInset, {
      target: document.body,
      props: {
        element: { ...mock_Fe, name: `` },
        elem_counts: { Fe: 10 },
      },
    })

    expect(document.querySelector(`strong`)?.textContent?.trim()).toBe(``)
  })
})
