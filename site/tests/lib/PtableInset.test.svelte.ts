import PtableInset from '$lib/PtableInset.svelte'
import type { ChemicalElement, ElementSymbol } from 'matterviz'
import { mount } from 'svelte'
import { describe, expect, it } from 'vitest'

type ElemCounts = Record<ElementSymbol, number>
const mock_Fe = { symbol: `Fe`, name: `Iron`, number: 26 } as ChemicalElement
const mock_H = { symbol: `H`, name: `Hydrogen`, number: 1 } as ChemicalElement

describe(`PtableInset.svelte`, () => {
  it(`renders element name and count from record`, () => {
    mount(PtableInset, {
      target: document.body,
      props: { element: mock_Fe, elem_counts: { Fe: 150, O: 50 } as ElemCounts },
    })

    const strong = document.querySelector(`strong`)
    expect(strong?.textContent).toContain(`Iron`)
    expect(strong?.textContent).toContain(`150`)
  })

  it(`renders element count from array using atomic number index`, () => {
    const counts = new Array(120).fill(0)
    counts[0] = 200 // H (number 1, index 0)

    mount(PtableInset, {
      target: document.body,
      props: { element: mock_H, elem_counts: counts },
    })

    const strong = document.querySelector(`strong`)
    expect(strong?.textContent).toContain(`Hydrogen`)
    expect(strong?.textContent).toContain(`200`)
  })

  it(`shows percentage by default and hides when show_percent=false`, () => {
    mount(PtableInset, {
      target: document.body,
      props: { element: mock_Fe, elem_counts: { Fe: 50, O: 50 } as ElemCounts },
    })
    expect(document.querySelector(`strong`)?.textContent).toContain(`%`)

    document.body.innerHTML = ``
    mount(PtableInset, {
      target: document.body,
      props: {
        element: mock_Fe,
        elem_counts: { Fe: 50, O: 50 } as ElemCounts,
        show_percent: false,
      },
    })
    expect(document.querySelector(`strong`)?.textContent).not.toContain(`%`)
  })

  it(`displays unit and renders HTML in unit prop`, () => {
    mount(PtableInset, {
      target: document.body,
      props: {
        element: mock_Fe,
        elem_counts: { Fe: 100 } as ElemCounts,
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
        elem_counts: { Fe: 100 } as ElemCounts,
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
        element: { symbol: `X`, number: 999 } as unknown as ChemicalElement,
        elem_counts: { X: 10 } as unknown as ElemCounts,
      },
    })

    expect(document.querySelector(`strong`)?.textContent?.trim()).toBe(``)
  })
})
