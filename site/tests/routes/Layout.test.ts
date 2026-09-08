import { goto } from '$app/navigation'
import { MODELS } from '$lib/models.svelte'
import Layout from '$routes/+layout.svelte'
import { createRawSnippet, tick } from 'svelte'
import { expect, it, vi } from 'vitest'
import { doc_query, mount, mount_with_url } from '../index'

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
