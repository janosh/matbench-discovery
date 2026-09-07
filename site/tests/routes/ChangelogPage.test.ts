import { load } from '$routes/changelog/+page.server'
import Page from '$routes/changelog/+page.svelte'
import { expect, it } from 'vitest'
import { mount } from '../index'

it(`renders changelog Markdown with release headings and linked code`, async () => {
  mount(Page, { target: document.body, props: { data: await load() } })

  expect(document.querySelector(`h1`)?.textContent).toBe(`Changelog`)
  const release_link = document.querySelector(`h2 a`)
  expect(release_link?.textContent).toBe(`v1.3.1`)
  expect(release_link?.getAttribute(`href`)).toBe(
    `https://github.com/janosh/matbench-discovery/compare/v1.3.0...v1.3.1`,
  )
  expect(document.querySelector(`li a code`)?.textContent).toBe(`#138`)
})
