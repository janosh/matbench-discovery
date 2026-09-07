import { goto } from '$app/navigation'
import { MODELS } from '$lib/models.svelte'
import Layout from '$routes/+layout.svelte'
import { tick } from 'svelte'
import { expect, it } from 'vitest'
import { mount } from '../index'

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
