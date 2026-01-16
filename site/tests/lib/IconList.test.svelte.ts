import IconList from '$lib/IconList.svelte'
import { mount } from 'svelte'
import { describe, expect, it } from 'vitest'

describe(`IconList.svelte`, () => {
  it.each([
    { icons: [], desc: `empty array` },
    { icons: undefined, desc: `undefined` },
  ])(`renders nothing when icons is $desc`, ({ icons }) => {
    mount(IconList, { target: document.body, props: { icons } })
    expect(document.querySelectorAll(`span, img`).length).toBe(0)
  })

  it(`renders Icon components for items with icon: prefix`, () => {
    mount(IconList, {
      target: document.body,
      props: {
        icons: [
          { id: `icon:GitHub`, name: `GitHub` },
          { id: `icon:Download`, name: `Download` },
        ],
      },
    })

    const spans = document.querySelectorAll(`span`)
    expect(spans.length).toBe(2)
    expect(spans[0].getAttribute(`title`)).toBe(`GitHub`)
    expect(spans[0].querySelector(`svg`)).toBeDefined()
  })

  it(`renders img elements for items with src`, () => {
    mount(IconList, {
      target: document.body,
      props: { icons: [{ src: `/logo.png`, name: `Company Logo` }] },
    })

    const img = document.querySelector(`img`)
    expect(img?.getAttribute(`src`)).toBe(`/logo.png`)
    expect(img?.getAttribute(`alt`)).toBe(`Company Logo logo`)
    expect(img?.getAttribute(`title`)).toBe(`Company Logo`)
    expect(img?.getAttribute(`style`)).toContain(`height: 1em`)
  })

  it(`forwards props to elements`, () => {
    mount(IconList, {
      target: document.body,
      props: {
        icons: [{ id: `icon:GitHub`, name: `GitHub` }],
        class: `custom-class`,
      },
    })

    expect(document.querySelector(`span`)?.classList.contains(`custom-class`)).toBe(true)
  })

  it(`skips items without icon: prefix and without src`, () => {
    mount(IconList, {
      target: document.body,
      props: {
        icons: [
          { id: `not-an-icon`, name: `Invalid` },
          { id: `icon:GitHub`, name: `GitHub` },
        ],
      },
    })

    expect(document.querySelectorAll(`span`).length).toBe(1)
  })
})
