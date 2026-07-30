import SelectToggle from '$lib/SelectToggle.svelte'
import { tick } from 'svelte'
import { describe, expect, it } from 'vitest'
import { mount } from '../index'

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

      // one button per option, only the selected one pressed
      const buttons = document.querySelectorAll(`button`)
      expect(
        [...buttons].map((button) => button.getAttribute(`aria-pressed`) === `true`),
      ).toStrictEqual(options.map((_, idx) => idx === expectedActiveIndex))

      // Check HTML rendering if expected
      expect(document.querySelector(`b`)?.textContent ?? null).toBe(
        expectHtml ? `Bold` : null,
      )
    },
  )

  it(`updates aria-pressed on the same buttons when one is clicked`, async () => {
    mount(SelectToggle, {
      target: document.body,
      props: {
        selected: `option1`,
        options: [
          { value: `option1`, label: `Option 1` },
          { value: `option2`, label: `Option 2` },
        ],
      },
    })

    // captured before the click so the assertions prove the existing nodes update
    // in place rather than the each block re-creating them
    const buttons = [...document.querySelectorAll(`button`)]
    const pressed = () => buttons.map((button) => button.getAttribute(`aria-pressed`))
    expect(pressed()).toStrictEqual([`true`, `false`])

    buttons[1].click()
    await tick()

    expect(pressed()).toStrictEqual([`false`, `true`])
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
      // interactive content inside a <button> is an invalid HTML content model
      expect(links[0]?.closest(`button`) ?? null).toBeNull()
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
