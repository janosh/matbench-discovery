import { afterNavigate } from '$app/navigation'
import { page } from '$app/state'
import { MODELS } from '$lib'
import {
  bind_comparison_url,
  compare_cells,
  comparison,
  mark_compared_rows,
  type CompareRow,
} from '$lib/model-comparison.svelte'
import ModelComparison from '$lib/model/ModelComparison.svelte'
import type { ModelData } from '$lib/types'
import { flushSync, tick } from 'svelte'
import { SvelteURL } from 'svelte/reactivity'
import { beforeEach, describe, expect, it, test, vi } from 'vitest'
import { get_scatter_plot_props, mount } from '../index'

// happy-dom never measures the plot, so capture ScatterPlot's props instead of its SVG
const plot_mocks = vi.hoisted(() => ({ ScatterPlot: vi.fn() }))
vi.mock(`matterviz`, async (import_original) => ({
  ...(await import_original<Record<string, unknown>>()),
  ScatterPlot: plot_mocks.ScatterPlot,
}))
type PlotProps = {
  series: { x: number[]; markers: string; point_style?: { radius: number }[] }[]
  x_axis: { label: string; scale_type: string }
  y_axis: { label: string }
}

const [key_a, key_b, key_c] = MODELS.map((model) => model.model_key)
// minimal stand-ins: a row whose accessor reads a single field off the "model"
const as_models = (vals: unknown[]) =>
  vals.map((val) => ({ val }) as unknown as ModelData)
const row = (better: `lower` | `higher` | null): CompareRow => ({
  key: `val`,
  label: `Val`,
  description: ``,
  better,
  value: (model) => (model as unknown as { val: unknown }).val,
})

beforeEach(() => {
  comparison.keys.clear()
  comparison.open = false
})

describe(`compare_cells`, () => {
  const field = as_models([0.5, 1, 2, 3, 7, null])
  test.each([
    // ranks numbers against the field, flags the best compared value, passes text through
    [
      `lower`,
      [3, 1, undefined, `text`, ``],
      [
        { text: `3`, rank: 4, n: 5, best: false },
        { text: `1`, rank: 2, n: 5, best: true },
        { text: `–` },
        { text: `text` },
        { text: `–` },
      ],
    ],
    // higher=better flips the rank; a lone numeric value is never "best"
    [`higher`, [3, `x`], [{ text: `3`, rank: 2, n: 5, best: false }, { text: `x` }]],
    // rows without a direction get neither rank nor best
    [null, [3, 1], [{ text: `3` }, { text: `1` }]],
  ] as const)(`better=%s -> %j`, (better, values, expected) => {
    expect(compare_cells(row(better), as_models([...values]), field)).toEqual(expected)
  })
})

describe(`comparison store`, () => {
  it(`toggles, replaces (dropping unknown keys, keeping order) and marks table rows`, () => {
    comparison.toggle(key_b)
    comparison.toggle(key_a)
    expect([...comparison.keys]).toEqual([key_b, key_a])
    expect(comparison.models.map((model) => model.model_key)).toEqual([key_b, key_a])
    comparison.toggle(key_b)
    expect([...comparison.keys]).toEqual([key_a])

    comparison.set([key_c, `no-such-model`, key_a, key_c])
    expect([...comparison.keys]).toEqual([key_c, key_a])

    const rows = () => [{ model_key: key_a }, { model_key: key_b }]
    expect(mark_compared_rows(rows(), false)).toEqual([
      { model_key: key_a, class: `highlight` },
      { model_key: key_b, class: undefined },
    ])
    expect(mark_compared_rows(rows(), true)).toEqual([
      { model_key: key_a, class: undefined },
    ])
  })
})

describe(`bind_comparison_url`, () => {
  // simulate SvelteKit: page.url is reactive there, so mutate one SvelteURL in place, mirror
  // it into location, then fire the afterNavigate hooks registered by bind_comparison_url
  const hooks = vi.mocked(afterNavigate).mock.calls
  const navigate = (url: string, type: `enter` | `link`, from_hook = 0) => {
    page.url.href = new URL(url, `http://localhost`).href
    history.replaceState(null, ``, `${page.url.pathname}${page.url.search}`)
    for (const [callback] of hooks.slice(from_hook)) callback({ type } as never)
    flushSync()
  }

  it(`reads shared links, keeps the selection across in-app navigation and re-applies it`, () => {
    Object.assign(page, { url: new SvelteURL(`http://localhost/`) }) // reactive like SvelteKit's
    const n_hooks = hooks.length
    const cleanup = $effect.root(() => bind_comparison_url())
    navigate(`/?compare=${key_a},no-such-model,${key_b}`, `enter`, n_hooks)
    expect([...comparison.keys]).toEqual([key_a, key_b])
    expect(comparison.open).toBe(true) // shared link with >1 model opens the drawer
    expect(location.search).toBe(`?compare=${key_a},${key_b}`) // stale key dropped

    navigate(`/tasks/md`, `link`, n_hooks) // absent param: selection follows the user, drawer closes
    expect([...comparison.keys]).toEqual([key_a, key_b])
    expect(comparison.open).toBe(false)
    expect(location.search).toBe(`?compare=${key_a},${key_b}`)

    comparison.toggle(key_b)
    flushSync()
    expect(location.search).toBe(`?compare=${key_a}`)
    comparison.keys.clear()
    flushSync()
    expect(location.search).toBe(``)

    navigate(`/?compare=${key_c}`, `enter`, n_hooks) // a single model is nothing to compare yet
    expect([...comparison.keys]).toEqual([key_c])
    expect(comparison.open).toBe(false)
    cleanup()
  })
})

describe(`ModelComparison drawer`, () => {
  it(`opens with a single selected model without looping and lists it`, async () => {
    comparison.toggle(key_a)
    mount(ModelComparison, { target: document.body })
    await tick()
    document.querySelector<HTMLButtonElement>(`.launcher button`)?.click()
    // a few flushes: a reactive write-back loop would hang here instead of settling
    for (let idx = 0; idx < 5; idx++) await tick()

    expect(comparison.open).toBe(true)
    expect([...comparison.keys]).toEqual([key_a])
    const header_cells = [...document.querySelectorAll(`thead th a`)]
    expect(header_cells.map((el) => el.textContent)).toEqual([MODELS[0].model_name])
    expect(document.querySelector(`.launcher`)).toBeNull() // hidden while the drawer is open
  })

  it(`removes models via the header buttons and clears via the launcher`, async () => {
    comparison.set([key_a, key_b, key_c])
    comparison.open = true
    mount(ModelComparison, { target: document.body })
    await tick()
    expect(document.querySelectorAll(`thead th a`)).toHaveLength(3)

    document
      .querySelector<HTMLButtonElement>(`thead th button[aria-label^="Remove"]`)
      ?.click()
    await tick()
    expect([...comparison.keys]).toEqual([key_b, key_c])
    expect(document.querySelectorAll(`thead th a`)).toHaveLength(2)

    comparison.open = false
    await tick()
    document
      .querySelector<HTMLButtonElement>(
        `.launcher button[aria-label="Clear model comparison"]`,
      )
      ?.click()
    await tick()
    expect(comparison.keys.size).toBe(0)
    expect(document.querySelector(`.launcher`)).toBeNull()
  })

  it(`picks scatter axes via selects, y defaulting to the current task page's metric`, async () => {
    Object.assign(page, { url: new URL(`http://localhost/tasks/phonons`) })
    const with_kappa = MODELS.filter((model) => model.metrics?.phonons?.kappa_103?.κ_SRME)
    comparison.set(with_kappa.slice(0, 2).map((model) => model.model_key))
    comparison.open = true
    mount(ModelComparison, { target: document.body })
    await tick()
    const [x_select, y_select] =
      document.querySelectorAll<HTMLSelectElement>(`dialog select`)
    const selected = (select: HTMLSelectElement) => select.selectedOptions[0]?.textContent
    expect(selected(x_select)).toBe(`Params`)
    expect(selected(y_select)).toBe(`κSRME`)
    const plot = () => get_scatter_plot_props(plot_mocks.ScatterPlot) as PlotProps
    expect(plot().x_axis).toMatchObject({ label: `Params`, scale_type: `log` })
    expect(plot().y_axis.label).toBe(`κ<sub>SRME</sub>`)
    // grey field + dashed Pareto frontier; compared models drawn last and larger
    const [field, frontier] = plot().series
    expect(frontier.markers).toBe(`line`)
    expect(field.point_style?.map((style) => style.radius).slice(-2)).toEqual([7, 7])

    const pick = (select: HTMLSelectElement, text: string) => {
      select.value =
        [...select.options].find((opt) => opt.textContent === text)?.value ?? ``
      select.dispatchEvent(new Event(`change`, { bubbles: true })) // Svelte delegates onchange
    }
    pick(x_select, `Training Materials`)
    pick(y_select, `F1`)
    await tick()
    expect(plot().x_axis.label).toBe(`Training Materials`)
    expect(plot().y_axis.label).toBe(`F1`)
  })
})
