import { load } from '#routes/api/+page.server.js'
import ApiPage from '#routes/api/+page.svelte'
import type * as node_fs from 'node:fs'
import { assert_ok, create_markdown } from 'svelte-widgets/markdown'
import { expect, it, vi } from 'vite-plus/test'
import { doc_query, mount_with_url } from '../index.js'

const env = vi.hoisted(() => ({ building: false }))
vi.mock(`$app/env`, () => ({
  get building() {
    return env.building
  },
}))
vi.mock(`node:fs`, async (import_original) => ({
  ...(await import_original<typeof node_fs>()),
  existsSync: () => false, // no generated API docs
}))

it.each([true, false])(
  `load without API docs throws only in builds (%s)`,
  async (building) => {
    env.building = building
    if (building) await expect(load()).rejects.toThrow(`No API docs`)
    else expect(await load()).toEqual({ modules: [] })
  },
)

it(`keeps API docs styles from leaking onto other pages' elements`, async () => {
  document.body.innerHTML = `<hr id="outside" /><h1 id="other"><kbd>key</kbd></h1>`
  const modules = [`<hr /><h1><kbd>key</kbd></h1>`]
  await mount_with_url(ApiPage, `http://localhost/api`, { props: { data: { modules } } })
  const style_of = (selector: string) => getComputedStyle(doc_query(selector))

  expect(style_of(`.api > hr`).marginBottom).toBe(`32px`) // 2em
  expect(style_of(`#outside`).marginBottom).not.toBe(`32px`)
  expect(style_of(`.api kbd`).fontSize).not.toBe(style_of(`#other > kbd`).fontSize)
})

it(`renders module docs as plain Markdown, so docstring braces stay literal`, async () => {
  const html = assert_ok(
    await create_markdown().render(`run_metadata ({} if absent) and \`{x}\``),
  )
  const rendered = document.createElement(`div`)
  rendered.innerHTML = html
  expect(rendered.textContent).toContain(`({} if absent)`)
  expect(rendered.querySelector(`code`)?.textContent).toBe(`{x}`)
})
