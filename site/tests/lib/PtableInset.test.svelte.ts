import PtableInset from '#lib/PtableInset.svelte'
import type { ChemicalElement } from 'matterviz'
import type { ComponentProps } from 'svelte'
import { describe, expect, it } from 'vite-plus/test'
import { doc_query, mount } from '../index.js'

const mock_Fe = { symbol: `Fe`, name: `Iron`, number: 26 } as ChemicalElement

describe(`PtableInset.svelte`, () => {
  it.each([
    [{ Fe: 150, O: 50 }, `Iron: 150 (75%)`],
    // elements missing from the counts show 0 rather than undefined
    [{ O: 50 }, `Iron: 0 (0%)`],
  ])(`renders element name, count and share from %o`, (elem_counts, expected) => {
    mount(PtableInset, {
      target: document.body,
      props: { element: mock_Fe, elem_counts },
    })
    expect(doc_query(`strong`).textContent?.replaceAll(/\s+/g, ` `).trim()).toBe(expected)
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

  it.each<{
    desc: string
    props: Partial<ComponentProps<typeof PtableInset>>
    read: (strong: HTMLElement) => unknown
    expected: unknown
  }>([
    {
      desc: `renders HTML in the unit prop`,
      props: { unit: `<sub>2</sub>`, show_percent: false },
      read: (strong) => strong.innerHTML.includes(`<sub>2</sub>`),
      expected: true,
    },
    {
      desc: `forwards class and style props`,
      props: { class: `custom-class`, style: `color: red;` },
      read: (strong) => [strong.className.includes(`custom-class`), strong.style.color],
      expected: [true, `red`],
    },
    {
      desc: `renders nothing when the element has no name`,
      props: { element: { ...mock_Fe, name: `` } },
      read: (strong) => strong.textContent?.trim(),
      expected: ``,
    },
  ])(`$desc`, ({ props, read, expected }) => {
    mount(PtableInset, {
      target: document.body,
      props: { element: mock_Fe, elem_counts: { Fe: 100 }, ...props },
    })
    expect(read(doc_query(`strong`))).toStrictEqual(expected)
  })
})
