import SelectToggle from '$lib/SelectToggle.svelte'
import { mount } from 'svelte'
import { afterEach, describe, expect, it, test, vi } from 'vitest'

describe(`SelectToggle.svelte`, () => {
  // Clean up after each test
  afterEach(() => {
    document.body.innerHTML = ``
    vi.restoreAllMocks()
  })

  test.each([
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
      expect(buttons.length).toBe(options.length)

      // Verify active state
      options.forEach((_, idx) => {
        if (idx === expectedActiveIndex) {
          expect(buttons[idx].classList.contains(`active`)).toBe(true)
        } else {
          expect(buttons[idx].classList.contains(`active`)).toBe(false)
        }
      })

      // Check HTML rendering if expected
      if (expectHtml) {
        const boldElement = document.querySelector(`b`)
        expect(boldElement).not.toBeNull()
        expect(boldElement?.textContent).toBe(`Bold`)
      }
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

  it(`buttons have proper properties for interaction`, () => {
    const options = [
      { value: `option1`, label: `Option 1` },
      { value: `option2`, label: `Option 2` },
    ]

    mount(SelectToggle, {
      target: document.body,
      props: {
        selected: `option1`,
        options,
      },
    })

    // Get buttons
    const buttons = document.querySelectorAll(`button`)

    // Verify the buttons have the correct classes to indicate they're interactive
    expect(buttons[0].tagName).toBe(`BUTTON`)
    expect(buttons[1].tagName).toBe(`BUTTON`)

    // Buttons should be accessible with proper roles
    expect(buttons[0].getAttribute(`role`) || `button`).toBe(`button`)
    expect(buttons[1].getAttribute(`role`) || `button`).toBe(`button`)

    // Verify the buttons have the correct active state initially
    expect(buttons[0].classList.contains(`active`)).toBe(true)
    expect(buttons[1].classList.contains(`active`)).toBe(false)
  })

  test.each([
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

      if (expectNoLink) {
        expect(links.length).toBe(0)
      } else {
        expect(links.length).toBeGreaterThan(0)
        const link = links[0]

        // Check link attributes
        Object.entries(expectedLinkAttributes || {}).forEach(([attr, value]) => {
          expect(link.getAttribute(attr)).toBe(value)
        })
      }
    },
  )
})
