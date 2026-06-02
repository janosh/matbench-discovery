import SelectToggle from '$lib/SelectToggle.svelte'
import { mount } from 'svelte'
import { describe, expect, it } from 'vitest'

describe(`SelectToggle.svelte`, () => {
  it.each([
    {
      name: `basic options`,
      options: [
        { value: `option1`, label: `Option 1`, tooltip: `First option` },
        { value: `option2`, label: `Option 2`, tooltip: `Second option` },
      ],
      selected: `option1`,
      expectedActiveIndex: 0,
    },
    {
      name: `with HTML content`,
      options: [
        { value: `option1`, label: `Option <b>Bold</b>` },
        { value: `option2`, label: `Option 2` },
      ],
      selected: `option1`,
      expectedActiveIndex: 0,
      expectHtml: true,
    },
    {
      name: `different selection`,
      options: [
        { value: `option1`, label: `Option 1` },
        { value: `option2`, label: `Option 2` },
      ],
      selected: `option2`,
      expectedActiveIndex: 1,
    },
  ])(
    `renders $name correctly`,
    ({ options, selected, expectedActiveIndex, expectHtml }) => {
      mount(SelectToggle, {
        target: document.body,
        props: { selected, options },
      })

      // Check that all options are rendered as buttons
      const buttons = document.querySelectorAll(`button`)
      expect(buttons).toHaveLength(options.length)

      // Verify active state
      expect(
        [...buttons].map((button) => button.classList.contains(`active`)),
      ).toStrictEqual(options.map((_, idx) => idx === expectedActiveIndex))

      // Check HTML rendering if expected
      expect(document.querySelector(`b`)?.textContent ?? null).toBe(
        expectHtml ? `Bold` : null,
      )
    },
  )

  it(`updates button active state when selection changes`, () => {
    const options = [
      { value: `option1`, label: `Option 1` },
      { value: `option2`, label: `Option 2` },
    ]

    // Start with option1 selected
    let selected = `option1`

    mount(SelectToggle, {
      target: document.body,
      props: { selected, options },
    })

    // Verify initial active state
    let buttons = document.querySelectorAll(`button`)
    expect(buttons[0].classList.contains(`active`)).toBe(true)
    expect(buttons[1].classList.contains(`active`)).toBe(false)

    // Remount with option2 selected to simulate state change
    document.body.innerHTML = ``
    selected = `option2`

    mount(SelectToggle, {
      target: document.body,
      props: { selected, options },
    })

    // Verify updated active state
    buttons = document.querySelectorAll(`button`)
    expect(buttons[0].classList.contains(`active`)).toBe(false)
    expect(buttons[1].classList.contains(`active`)).toBe(true)
  })

  it.each([
    {
      name: `with link`,
      options: [{ value: `option1`, label: `Option 1`, link: `https://example.com` }],
      selected: `option1`,
      expectedLinkAttributes: {
        href: `https://example.com`,
        target: `_blank`,
        rel: `noopener noreferrer`,
      },
    },
    {
      name: `without link`,
      options: [{ value: `option1`, label: `Option 1` }],
      selected: `option1`,
      expectNoLink: true,
    },
  ])(
    `renders $name correctly`,
    ({ options, selected, expectedLinkAttributes, expectNoLink }) => {
      mount(SelectToggle, {
        target: document.body,
        props: { selected, options },
      })

      const links = document.querySelectorAll(`a`)

      expect(links).toHaveLength(expectNoLink ? 0 : 1)
      expect(links[0]?.getAttribute(`href`) ?? null).toBe(
        expectedLinkAttributes?.href ?? null,
      )
      expect(links[0]?.getAttribute(`target`) ?? null).toBe(
        expectedLinkAttributes?.target ?? null,
      )
      expect(links[0]?.getAttribute(`rel`) ?? null).toBe(
        expectedLinkAttributes?.rel ?? null,
      )
    },
  )
})
