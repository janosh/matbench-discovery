import { load } from '#routes/changelog/+page.server.js'
import Page from '#routes/changelog/+page.svelte'
import { expect, it } from 'vite-plus/test'
import { doc_query, mount } from '../index.js'

it(`renders changelog Markdown with release headings and linked code`, async () => {
  mount(Page, { target: document.body, props: { data: await load() } })

  expect(document.querySelector(`h1`)?.textContent).toBe(`Changelog`)
  // newest release first; match the shape, not a pinned version, so releases don't break it
  const release_link = doc_query<HTMLAnchorElement>(`h2 a`)
  const version = release_link.textContent ?? ``
  expect(version).toMatch(/^v\d+\.\d+\.\d+$/)
  expect(release_link.getAttribute(`href`)).toMatch(
    new RegExp(
      `^https://github\\.com/janosh/matbench-discovery/compare/v[\\d.]+\\.\\.\\.${version}$`,
    ),
  )
  expect(doc_query(`li a code`).textContent).toMatch(/^#\d+$/)
})
