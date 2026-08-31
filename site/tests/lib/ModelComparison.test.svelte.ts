import { afterNavigate } from '$app/navigation'
import { page } from '$app/state'
import { DATASETS, MODELS } from '$lib'
import {
  bind_comparison_url,
  COMPARE_GROUPS,
  compare_cells,
  comparison,
  mark_compared_rows,
  type CompareRow,
} from '$lib/model-comparison.svelte'
import CompareToggle from '$lib/model/CompareToggle.svelte'
import ModelComparison from '$lib/model/ModelComparison.svelte'
import type { ModelData } from '$lib/types'
import { flushSync, tick } from 'svelte'
import { afterEach, beforeEach, describe, expect, it, test, vi } from 'vitest'
import { doc_query, get_scatter_plot_props, mount } from '../index'

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
  // compared models are field members where the value exists there (like active models on
  // the real leaderboard), else stand-ins outside the field (superseded/deprecated models)
  const pick = (vals: readonly unknown[]) =>
    vals.map(
      (val) =>
        field.find((model) => (model as unknown as { val: unknown }).val === val) ??
        as_models([val])[0],
    )
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
    // a model outside the field is ranked against the field plus itself (never rank > n)
    [
      `lower`,
      [10, 0.1],
      [
        { text: `10`, rank: 6, n: 6, best: false },
        { text: `0.1`, rank: 1, n: 6, best: true },
      ],
    ],
    // rows without a direction get neither rank nor best
    [null, [3, 1], [{ text: `3` }, { text: `1` }]],
  ] as const)(`better=%s -> %j`, (better, values, expected) => {
    expect(compare_cells(row(better), pick(values), field)).toEqual(expected)
  })

  it(`renders text rows as parts and attaches titles to numeric cells`, () => {
    const parts_row: CompareRow = {
      ...row(null),
      parts: (model) => {
        const val = (model as unknown as { val: string }).val
        return val ? [{ text: val, href: `/data/${val}` }, { text: ` + ` }] : []
      },
    }
    expect(compare_cells(parts_row, as_models([`mptrj`, ``, 3]), field)).toEqual([
      {
        text: `mptrj + `,
        parts: [{ text: `mptrj`, href: `/data/mptrj` }, { text: ` + ` }],
      },
      { text: `–` }, // empty parts hide the cell rather than rendering an empty link
      { text: `3` }, // numbers still format even on a row with `parts`
    ])
    const titled_row: CompareRow = {
      ...row(`lower`),
      title: (model) =>
        (model as unknown as { val: number }).val > 1 ? `big` : undefined,
    }
    expect(compare_cells(titled_row, pick([3, 1]), field)).toEqual([
      { text: `3`, title: `big`, rank: 4, n: 5, best: false },
      { text: `1`, rank: 2, n: 5, best: true },
    ])
  })

  it(`links training sets, papers, repos and authors of real models`, () => {
    const model = MODELS.find(
      (md) => md.paper && md.repo && md.authors[0].url && md.authors.length > 3,
    )
    if (!model) throw new Error(`no model with paper, repo, author URL and >3 authors`)
    const rows = Object.fromEntries(
      COMPARE_GROUPS[0].rows.map((cmp_row) => [
        cmp_row.key,
        compare_cells(cmp_row, [model])[0],
      ]),
    )
    expect(rows.training_sets.parts?.filter((part) => part.href)).toEqual(
      model.training_sets.map((key) => ({
        text: key,
        href: `/data/${DATASETS[key].slug}`,
        title: expect.stringContaining(`${DATASETS[key].name}: `),
      })),
    )
    const links = Object.fromEntries(
      rows.links.parts
        ?.filter((part) => part.href)
        .map((part) => [part.text, part.href]) ?? [],
    )
    expect(links).toMatchObject({ Paper: model.paper, Repo: model.repo })
    expect(rows.authors.parts?.[0]).toEqual({
      text: model.authors[0].name,
      href: model.authors[0].url,
      title: model.authors[0].affiliation,
    })
    expect(rows.authors.parts?.at(-1)?.text).toBe(`+${model.authors.length - 3} more`)
    expect(rows.targets.parts).toEqual([
      { text: model.targets, title: expect.stringMatching(/^Energy/) },
    ])
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
  // simulate SvelteKit: its page.url is a $state.raw holding a fresh URL object per navigation
  // (kit/src/runtime/client/client.js), so swap the object, mirror it into location, then
  // fire the afterNavigate hooks registered by bind_comparison_url
  let kit_url = $state.raw(new URL(`http://localhost/`))
  const hooks = vi.mocked(afterNavigate).mock.calls
  const navigate = (url: string, type: `enter` | `link`, from_hook = 0) => {
    kit_url = new URL(url, `http://localhost`)
    history.replaceState(null, ``, `${kit_url.pathname}${kit_url.search}`)
    for (const [callback] of hooks.slice(from_hook)) callback({ type } as never)
    flushSync()
  }
  let cleanup = () => {}

  beforeEach(() => {
    Object.defineProperty(page, `url`, { get: () => kit_url, configurable: true })
  })
  afterEach(() => {
    cleanup() // tear down the effect root even when an assertion above it failed
    // back to the plain, assignable mock property other tests (and tests/index) rely on
    Object.defineProperty(page, `url`, {
      value: new URL(`http://localhost/`),
      writable: true,
      enumerable: true,
      configurable: true,
    })
  })

  it(`reads shared links, keeps the selection across in-app navigation and re-applies it`, () => {
    const n_hooks = hooks.length
    cleanup = $effect.root(() => bind_comparison_url())
    navigate(`/?compare=${key_a},no-such-model,${key_b}`, `enter`, n_hooks)
    expect([...comparison.keys]).toEqual([key_a, key_b])
    expect(comparison.open).toBe(true) // shared link with >1 model opens the dialog
    expect(location.search).toBe(`?compare=${key_a},${key_b}`) // stale key dropped

    navigate(`/tasks/md`, `link`, n_hooks) // absent param: selection follows the user, dialog closes
    expect([...comparison.keys]).toEqual([key_a, key_b])
    expect(comparison.open).toBe(false)
    expect(location.search).toBe(`?compare=${key_a},${key_b}`)

    // same-path navigation that only drops the query (nav link of the current page) must
    // re-apply the param too
    navigate(`/tasks/md`, `link`, n_hooks)
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
  })
})

describe(`ModelComparison dialog`, () => {
  it(`opens with a single selected model without looping and lists it`, async () => {
    comparison.toggle(key_a)
    mount(ModelComparison, { target: document.body })
    await tick()
    // a lone pick shows the tray with its chip and a hint about the next step
    const tray_text = document.querySelector(`.tray`)?.textContent ?? ``
    expect(tray_text).toContain(MODELS[0].model_name)
    expect(tray_text).toContain(`Pick a 2nd model`)
    expect(tray_text).toContain(`Open comparison`)
    document.querySelector<HTMLButtonElement>(`.tray button.open`)?.click()
    // a few flushes: a reactive write-back loop would throw effect_update_depth_exceeded
    // here instead of settling
    for (let idx = 0; idx < 5; idx++) await tick()

    expect(comparison.open).toBe(true)
    expect([...comparison.keys]).toEqual([key_a])
    const header_cells = [...document.querySelectorAll(`thead th a`)]
    expect(header_cells.map((el) => el.textContent)).toEqual([MODELS[0].model_name])
    expect(document.querySelector(`.tray`)).toBeNull() // hidden while the dialog is open
    // with nothing to compare yet, the model picker is the next step and takes focus
    expect(document.activeElement).toBe(
      document.querySelector(`dialog input[role="combobox"]`),
    )
  })

  it(`tray lists chips with remove buttons and counts the rest`, async () => {
    const keys = MODELS.slice(0, 6).map((model) => model.model_key)
    comparison.set(keys)
    mount(ModelComparison, { target: document.body })
    await tick()
    const chips = () => [...document.querySelectorAll(`.tray li`)]
    expect(chips().map((li) => li.textContent?.trim())).toEqual([
      ...MODELS.slice(0, 4).map((model) => `● ${model.model_name}`),
      `+2 more`,
    ])
    expect(document.querySelector(`.tray button.open`)?.textContent).toContain(
      `Compare 6 models`,
    )
    expect(document.querySelector(`.tray .hint`)).toBeNull()
    chips()[1].querySelector(`button`)?.click()
    await tick()
    expect([...comparison.keys]).toEqual(keys.filter((key) => key !== keys[1]))
    expect(chips()).toHaveLength(5) // 4 chips + "+1 more"
  })

  it(`CompareToggle in open_dialog mode adds the model and opens the dialog`, async () => {
    mount(CompareToggle, {
      target: document.body,
      props: { model_key: key_a, open_dialog: true },
    })
    const button = doc_query<HTMLButtonElement>(`button`)
    expect(button.textContent?.trim()).toBe(`Compare with…`)
    button.click()
    await tick()
    expect([...comparison.keys]).toEqual([key_a])
    expect(comparison.open).toBe(true)
    // a second click keeps the model selected rather than toggling it away
    comparison.toggle(key_b)
    await tick()
    expect(button.textContent?.trim()).toBe(`Comparing with 1 other`)
    button.click()
    await tick()
    expect([...comparison.keys]).toEqual([key_a, key_b])
  })

  it(`removes models via the header buttons and clears via the tray`, async () => {
    comparison.set([key_a, key_b, key_c])
    comparison.open = true
    mount(ModelComparison, { target: document.body })
    await tick()
    expect(document.querySelectorAll(`thead th a`)).toHaveLength(3)
    // dataset links stay in-app, paper/repo open new tabs
    const data_links =
      document.querySelectorAll<HTMLAnchorElement>(`td a[href^="/data/"]`)
    expect(data_links.length).toBeGreaterThan(0)
    expect(data_links[0].target).toBe(``)
    const ext_link = document.querySelector<HTMLAnchorElement>(`td a[href^="http"]`)
    expect(ext_link?.target).toBe(`_blank`)

    const remove_first = () =>
      document
        .querySelector<HTMLButtonElement>(`thead th button[aria-label^="Remove"]`)
        ?.click()
    remove_first()
    await tick()
    expect([...comparison.keys]).toEqual([key_b, key_c])
    expect(document.querySelectorAll(`thead th a`)).toHaveLength(2)
    // dropping to a single model mid-session must not yank focus to the picker (that
    // only happens when the dialog is opened with <=1 models)
    remove_first()
    for (let idx = 0; idx < 3; idx++) await tick()
    expect([...comparison.keys]).toEqual([key_c])
    expect(document.activeElement).not.toBe(
      document.querySelector(`dialog input[role="combobox"]`),
    )

    // dropdown lists models by CPS (best first, gray score after the name) until the sort
    // toggle switches it to newest first
    const option_labels = () =>
      [...document.querySelectorAll(`dialog ul.options li`)].map((li) =>
        li.textContent?.replaceAll(/\s+/g, ` `).trim(),
      )
    const unselected = MODELS.filter((model) => !comparison.keys.has(model.model_key))
    const cps_of = (model: ModelData) => model.CPS ?? Number.NaN
    const by_cps = unselected
      .filter((model) => cps_of(model) > 0)
      .toSorted((md_1, md_2) => cps_of(md_2) - cps_of(md_1))
    expect(option_labels().slice(0, 2)).toEqual(
      by_cps
        .slice(0, 2)
        .map((model) => `${model.model_name} ${cps_of(model).toFixed(3)}`),
    )
    const sort_select = document.querySelector<HTMLSelectElement>(`.option-sort select`)
    if (!sort_select) throw new Error(`sort select not rendered`)
    sort_select.value = `added`
    sort_select.dispatchEvent(new Event(`change`, { bubbles: true }))
    await tick()
    const newest = unselected.toSorted((md_1, md_2) =>
      (md_2.dates.benchmark_added ?? ``).localeCompare(md_1.dates.benchmark_added ?? ``),
    )[0]
    expect(option_labels()[0]).toContain(`${newest.model_name} `)
    expect(option_labels()[0]).toContain(`· ${newest.dates.benchmark_added}`)

    // picking an option adds that model (MultiSelect's bind:selected writes back to the
    // store), and the dialog's own close button closes it
    document.querySelector<HTMLLIElement>(`dialog ul.options li`)?.click()
    await tick()
    expect([...comparison.keys]).toEqual([key_c, newest.model_key])
    expect(document.querySelectorAll(`thead th a`)).toHaveLength(2)
    document
      .querySelector<HTMLButtonElement>(`button[aria-label="Close model comparison"]`)
      ?.click()
    await tick()
    expect(comparison.open).toBe(false)
    document
      .querySelector<HTMLButtonElement>(
        `.tray button[aria-label="Clear model comparison"]`,
      )
      ?.click()
    await tick()
    expect(comparison.keys.size).toBe(0)
    expect(document.querySelector(`.tray`)).toBeNull()
  })

  it(`picks scatter axes via selects, y defaulting to the current task page's metric`, async () => {
    Object.assign(page, { url: new URL(`http://localhost/tasks/phonons`) })
    const with_kappa = MODELS.filter((model) => model.metrics?.phonons?.kappa_103?.κ_SRME)
    comparison.set(with_kappa.slice(0, 2).map((model) => model.model_key))
    comparison.open = true
    mount(ModelComparison, { target: document.body })
    await tick()
    const [x_select, y_select] = document.querySelectorAll<HTMLSelectElement>(
      `dialog :not(.option-sort) > select`,
    )
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
