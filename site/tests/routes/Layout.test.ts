import { goto } from '$app/navigation'
import { MODELS } from '$lib/models.svelte'
import { comparison } from '$lib/model-comparison.svelte'
import Layout from '$routes/+layout.svelte'
import { createRawSnippet, tick } from 'svelte'
import { expect, it, vi } from 'vitest'
import { doc_query, mount, mount_with_url } from '../index'

it(`loads comparison on demand and retains its controls between openings`, async () => {
  comparison.keys.clear()
  comparison.open = false
  await mount_with_url(Layout, `http://localhost/`)
  expect(document.querySelector(`dialog[aria-label="Model comparison"]`)).toBeNull()
  comparison.open_with(MODELS[0].model_key)
  const dialog = await vi.waitFor(
    () => doc_query<HTMLDialogElement>(`dialog[aria-label="Model comparison"]`),
    { timeout: 5000 },
  )
  const select = doc_query<HTMLSelectElement>(`select`, dialog)
  select.value = `added`
  select.dispatchEvent(new Event(`change`, { bubbles: true }))
  comparison.open = false
  await tick()
  comparison.open = true
  await tick()
  expect(
    doc_query<HTMLSelectElement>(`dialog[aria-label="Model comparison"] select`).value,
  ).toBe(`added`)
  comparison.open = false
  comparison.keys.clear()
})

it.each([`/api`, `/data`])(
  `shows the table of contents for three headings on %s`,
  async (route) => {
    await mount_with_url(Layout, `http://localhost${route}`, {
      props: {
        children: createRawSnippet(() => ({
          render: () =>
            `<section><h1>Overview</h1><h2>First</h2><h2>Second</h2></section>`,
        })),
      },
    })
    await vi.waitFor(() => {
      expect(document.querySelector(`.toc button`) !== null).toBe(route === `/api`)
    })
    if (route === `/data`) return

    expect(doc_query(`.toc`).style.zIndex).toBe(`1`)
    doc_query<HTMLButtonElement>(`.toc button`).click()
    await tick()
    expect(
      [...document.querySelectorAll(`.toc nav a`)].map((link) => link.textContent),
    ).toEqual([`Overview`, `First`, `Second`])
    expect(doc_query(`.toc nav`).style.fontSize).toBe(`7pt`)
    expect(doc_query(`.toc-title`).style.margin).toBe(`3pt`)
  },
)

it.each([`/data`, `/models/${MODELS[0].model_key}`])(
  `navigates to %s from the command menu`,
  async (route) => {
    mount(Layout, { target: document.body })
    await tick()
    window.dispatchEvent(new KeyboardEvent(`keydown`, { key: `k`, ctrlKey: true }))
    await tick()
    const menu = document.querySelector(`dialog[aria-label="Command menu"]`)
    const option = [...(menu?.querySelectorAll(`li[role="option"]`) ?? [])].find(
      (element) => element.textContent?.trim() === route,
    )
    expect(option).toBeDefined()
    option?.dispatchEvent(new MouseEvent(`click`, { bubbles: true }))
    await tick()
    expect(goto).toHaveBeenLastCalledWith(route)
  },
)
