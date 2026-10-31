import ApiPage from '#routes/api/+page.svelte'
import { assert_ok, create_markdown } from 'svelte-widgets/markdown'
import { expect, it } from 'vite-plus/test'
import { doc_query, mount_with_url } from '../index.js'

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
