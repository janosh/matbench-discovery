import { page } from '$app/state'
import { Nav } from '$lib'
import { mount } from 'svelte'
import { beforeEach, describe, expect, test, vi } from 'vitest'

describe(`Nav`, () => {
  beforeEach(() => {
    document.body.innerHTML = ``
  })

  describe(`rendering routes`, () => {
    test.each([
      {
        name: `with basic routes`,
        routes: [`/home`, `/about`, `/contact`],
        expected_hrefs: [`/home`, `/about`, `/contact`],
        expected_texts: [`/home`, `/about`, `/contact`],
      },
      {
        name: `with title/href pairs`,
        routes: [
          [`Home`, `/home`],
          [`About Us`, `/about`],
          [`Contact`, `/contact`],
        ] as [string, string][],
        expected_hrefs: [`/home`, `/about`, `/contact`],
        expected_texts: [`Home`, `About Us`, `Contact`],
      },
      {
        name: `with mixed string and pair routes`,
        routes: [`/home`, [`About Us`, `/about`] as [string, string], `/contact`],
        expected_hrefs: [`/home`, `/about`, `/contact`],
        expected_texts: [`/home`, `About Us`, `/contact`],
      },
      { name: `with empty array`, routes: [], expected_hrefs: [], expected_texts: [] },
    ])(`renders $name`, ({ routes, expected_hrefs, expected_texts }) => {
      mount(Nav, {
        target: document.body,
        props: { routes },
      })

      const links = document.querySelectorAll(`a`)
      expect(links.length).toBe(expected_hrefs.length)

      // Check links have correct hrefs and text
      expected_hrefs.forEach((href, idx) => {
        expect(links[idx]?.getAttribute(`href`)).toBe(href)
        expect(links[idx]?.textContent).toBe(expected_texts[idx])
      })
    })
  })

  describe(`separator bullets`, () => {
    test.each([
      { routes: [`/home`, `/about`, `/contact`], expected_bullets: 2 },
      { routes: [`/home`], expected_bullets: 0 },
      { routes: [`/home`, `/about`], expected_bullets: 1 },
      { routes: [], expected_bullets: 0 },
    ])(
      `applies $expected_bullets bullets for $routes.length routes`,
      ({ routes, expected_bullets }) => {
        mount(Nav, {
          target: document.body,
          props: { routes },
        })

        const bullets = document.querySelectorAll(`strong`)
        expect(bullets.length).toBe(expected_bullets)

        if (expected_bullets > 0) {
          // Check all bullets contain the bullet character
          Array.from(bullets).forEach((bullet) => {
            expect(bullet.textContent).toBe(`•`)
          })
        }
      },
    )
  })

  describe(`aria-current attribute`, () => {
    test.each([
      {
        path: `/about`,
        routes: [`/home`, `/about`, `/contact`],
        expected_aria_current: [null, `page`, null],
        desc: `exact match`,
      },
      {
        path: `/about/team`,
        routes: [`/home`, `/about`, `/contact`],
        expected_aria_current: [null, `page`, null],
        desc: `partial match`,
      },
      {
        path: `/contact/form`,
        routes: [`/home`, `/about`, `/contact`],
        expected_aria_current: [null, null, `page`],
        desc: `different partial match`,
      },
      {
        path: `/`,
        routes: [`/home`, `/about`, `/`],
        expected_aria_current: [null, null, `page`],
        desc: `root path exact match`,
      },
      {
        path: `/blog`,
        routes: [`/home`, `/about`, `/contact`],
        expected_aria_current: [null, null, null],
        desc: `no match`,
      },
    ])(
      `applies aria-current correctly for $desc`,
      ({ path, routes, expected_aria_current }) => {
        // Override the mock for this specific test
        vi.mocked(page).url.pathname = path

        mount(Nav, {
          target: document.body,
          props: { routes },
        })

        const links = document.querySelectorAll(`a`)

        for (let idx = 0; idx < expected_aria_current.length; idx++) {
          expect(links[idx]?.getAttribute(`aria-current`)).toBe(
            expected_aria_current[idx],
          )
        }
      },
    )
  })

  describe(`custom styling`, () => {
    test.each([
      {
        style: `color: red; font-size: 20px;`,
        expected_style: `color: red; font-size: 20px;`,
      },
      { style: null, expected_style: null },
      { style: ``, expected_style: `` },
    ])(`applies style correctly: $style`, ({ style, expected_style }) => {
      mount(Nav, {
        target: document.body,
        props: { routes: [`/home`, `/about`], style },
      })

      const nav = document.querySelector(`nav`)
      expect(nav?.getAttribute(`style`)).toBe(expected_style)
    })
  })
})
