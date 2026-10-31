import { afterNavigate, onNavigate, type AfterNavigate } from '$app/navigation'
import { page } from '$app/state'
import DATASETS from '$data/datasets.yml'
import { ALL_METRICS } from '#lib/labels.js'
import {
  ACTIVE_MODELS,
  bind_score_weights,
  MODELS,
  score_weight_records,
} from '#lib/models.svelte.js'
import {
  bind_comparison_url,
  COMPARE_GROUPS,
  compare_cells,
  comparison,
  format_duration,
  mark_compared_rows,
  type CompareRow,
} from '#lib/model-comparison.svelte.js'
import { bind_url_params } from '#lib/url-state.svelte.js'
import ModelComparison from '#lib/model/ModelComparison.svelte'
import ModelPage from '#routes/models/[slug]/+page.svelte'
import type { ModelData } from '#lib/types.js'
import { flushSync, tick } from 'svelte'
import {
  afterEach,
  beforeEach,
  describe,
  expect,
  it,
  onTestFinished,
  test,
  vi,
} from 'vite-plus/test'
import {
  doc_query,
  get_scatter_plot_props,
  mount,
  query_param,
  reset_comparison,
} from '../index.js'

// happy-dom never measures the plot, so capture ScatterPlot's props instead of its SVG
const plot_mocks = vi.hoisted(() => ({ ScatterPlot: vi.fn() }))
vi.mock(`matterviz/plot`, async (import_original) => ({
  ...(await import_original<Record<string, unknown>>()),
  ScatterPlot: plot_mocks.ScatterPlot,
}))
type PlotProps = {
  series: {
    id?: string
    markers: string
    point_style?: { stroke?: string; fill_opacity?: number }
    point_label?: { text: string }
  }[]
  x_axis: {
    label: string
    scale_type: string
  }
  y_axis: { label: string; scale_type: string }
  on_axis_change: (axis: `x` | `y`, key: string) => void
  point_events: {
    onclick: (payload: { point: { metadata: { model_key: string } } }) => void
  }
}

const [key_a, key_b, key_c] = MODELS.map((model) => model.model_key)
// minimal stand-ins: a row whose accessor reads a single field off the "model"
type StandIn = { val: unknown }
const as_models = (vals: unknown[]) =>
  vals.map((val) => ({ val }) as unknown as ModelData)
const val_of = (model: ModelData): unknown => (model as unknown as StandIn).val
const row = (better: `lower` | `higher` | null): CompareRow => ({
  key: `val`,
  label: `Val`,
  description: ``,
  better,
  value: val_of,
})
// fire the most recently registered afterNavigate callback
const last_after_navigate = (navigation: AfterNavigate): void =>
  vi.mocked(afterNavigate).mock.calls.at(-1)?.[0](navigation)

// SvelteKit's page.url is a $state.raw holding a fresh URL object per navigation
// (kit/src/runtime/client/client.js). Swap the plain mock property for a reactive
// accessor in this describe block; assigning page.url then re-runs dependent effects.
function with_reactive_page_url(): void {
  let kit_url = $state.raw<typeof page.url>(new URL(`http://localhost/`))
  beforeEach(() => {
    kit_url = page.url
    Object.defineProperty(page, `url`, {
      get: () => kit_url,
      set: (url: URL) => {
        kit_url = url
      },
      configurable: true,
    })
  })
  afterEach(() => {
    // back to the plain, assignable mock property other tests (and tests/index) rely on
    Object.defineProperty(page, `url`, {
      value: new URL(`http://localhost/`),
      writable: true,
      enumerable: true,
      configurable: true,
    })
  })
}

beforeEach(reset_comparison)

describe(`compare_cells`, () => {
  it.each([
    [0.4, `0s`],
    [52, `52s`],
    [108.15, `1m 48s`],
    [10_720.95, `2h 58m 41s`],
    [7_205, `2h 0m 5s`],
    [220_585.84, `2d 13h 16m`],
  ])(`formats %s s as %s`, (seconds, expected) => {
    expect(format_duration(seconds)).toBe(expected)
  })

  it(`shows units on values: durations for seconds, plain decimals for small GB`, () => {
    const [model] = MODELS.filter((entry) => entry.metrics?.md?.cost?.max_gpu_mem_gb)
    if (!model) throw new Error(`no model with MD VRAM`)
    const row_for = (key: string) =>
      COMPARE_GROUPS.flatMap((group) => group.rows).find((row) => row.key === key)
    const vram_row = row_for(ALL_METRICS.md_max_gpu_mem_gb.key)
    const time_row = row_for(ALL_METRICS.md_run_time_sec.key)
    if (!vram_row || !time_row) throw new Error(`missing cost rows`)
    expect(compare_cells(vram_row, [model])[0].text).toMatch(/^\d+(\.\d+)? GB$/)
    expect(compare_cells(time_row, [model])[0].text).toMatch(/^(\d+[dhms] ?){1,3}$/)
  })

  const field = as_models([0.5, 1, 2, 3, 7, null])
  // compared models are field members where the value exists there (like active models on
  // the real leaderboard), else stand-ins outside the field (superseded/deprecated models)
  const pick = (vals: readonly unknown[]) =>
    vals.map((val) => field.find((model) => val_of(model) === val) ?? as_models([val])[0])
  test.each([
    // ranks numbers against the field, flags the best compared value, passes text through
    [
      `lower`,
      [3, 1, undefined, `text`, ``],
      [
        { text: `3`, rank: 4, n: 5, best: false },
        { text: `1`, rank: 2, n: 5, best: true },
        { text: `–`, title: `Val: not reported.` },
        { text: `text` },
        { text: `–`, title: `Val: not reported.` },
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
        const val = val_of(model)
        return typeof val === `string` && val
          ? [{ text: val, href: `/data/${val}` }, { text: ` + ` }]
          : []
      },
    }
    expect(compare_cells(parts_row, as_models([`mptrj`, ``, 3]), field)).toEqual([
      {
        text: `mptrj + `,
        parts: [{ text: `mptrj`, href: `/data/mptrj` }, { text: ` + ` }],
      },
      { text: `–`, title: `Val: not reported.` }, // empty parts don't create empty links
      { text: `3` }, // numbers still format even on a row with `parts`
    ])
    const titled_row: CompareRow = {
      ...row(`lower`),
      title: (model) => (Number(val_of(model)) > 1 ? `big` : undefined),
    }
    expect(compare_cells(titled_row, pick([3, 1]), field)).toEqual([
      { text: `3`, title: `big`, rank: 4, n: 5, best: false },
      { text: `1`, rank: 2, n: 5, best: true },
    ])
    const gnome = MODELS.find(({ model_key }) => model_key === `gnome`)
    if (!gnome) throw new Error(`Missing GNoME test fixture`)
    expect(compare_cells(ALL_METRICS.RMSD, [gnome])[0]).toEqual({
      text: `–`,
      title: `Geometry Optimization: no results reported. Model weights are not publicly available.`,
    })
    expect(
      compare_cells(ALL_METRICS.RMSD, [
        {
          ...gnome,
          metrics: {
            geo_opt: {
              status: `not_available`,
              reason: `Missing <reference> & metadata.`,
            },
          },
        },
      ])[0].title,
    ).toContain(`Missing &lt;reference&gt; &amp; metadata.`)
  })

  it(`names the GPU and source of speed rows on hover`, () => {
    const with_md = (model_key: string, md: object) =>
      ({ ...MODELS[0], model_key, metrics: { md } }) as ModelData
    const models = [
      with_md(`authors`, {
        cost: { run_time_sec: 3e4, hardware: `A100 <80GB>`, source: `author_reported` },
        combined_score: 0.7,
      }),
      with_md(`untimed`, { cost: { hardware: `H200` } }),
    ]
    const all_rows = COMPARE_GROUPS.flatMap(({ rows }) => rows)
    const cell = (key: string) => {
      const cmp_row = all_rows.find((candidate) => candidate.key === key)
      if (!cmp_row) throw new Error(`No ${key} row`)
      return compare_cells(cmp_row, models, models)
    }
    // numeric-cell titles render as HTML, so YAML hardware strings arrive escaped
    expect(cell(`md_run_time_sec`)[0].title).toBe(
      `Timed on A100 &lt;80GB&gt; (author-reported), not normalized across hardware`,
    )
    expect(cell(`combined_score`)[0]).toMatchObject({
      text: `0.700`,
      title: expect.stringMatching(/^Speed subscore timed on A100 &lt;80GB&gt;/),
    })
    expect(cell(`md_hardware`).map(({ text }) => text)).toEqual([
      `A100 <80GB> (author-reported)`,
      `–`,
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

it(`reads shallow edits on late mount but the original router query on navigation`, () => {
  const navigation = { type: `link` } as AfterNavigate
  onTestFinished(
    $effect.root(() => {
      bind_url_params(null, () => [[`columns`, `Model,CPS`]])
    }),
  )
  last_after_navigate(navigation)
  flushSync()
  expect(location.search).toBe(`?columns=Model,CPS`)
  expect(page.url.search).toBe(``)

  const reads: (string | null)[] = []
  onTestFinished(
    $effect.root(() => {
      bind_url_params(
        (params) => {
          reads.push(params.get(`columns`))
        },
        () => [],
      )
    }),
  )
  flushSync()
  expect(reads).toEqual([`Model,CPS`])

  // Another binding has already written its automatic defaults. The navigation
  // callback must not mistake those writes for explicit columns from a shared URL.
  last_after_navigate(navigation)
  flushSync()
  expect(reads).toEqual([`Model,CPS`, null])

  // Back restores shallow edits in location, but Kit's router URL is still original.
  const popped = { type: `popstate`, complete: Promise.resolve() } as AfterNavigate
  history.replaceState(null, ``, `/?columns=Model,F1`)
  for (const [callback] of vi.mocked(onNavigate).mock.calls) {
    void callback(popped as Parameters<Parameters<typeof onNavigate>[0]>[0])
  }
  // A destination binding may write defaults before the remaining readers run.
  history.replaceState(null, ``, `/?columns=Model,CPS`)
  last_after_navigate(popped)
  flushSync()
  expect(reads.at(-1)).toBe(`Model,F1`)
  last_after_navigate(navigation)
  expect(reads.at(-1)).toBeNull()
})

describe(`bind_comparison_url`, () => {
  with_reactive_page_url()
  // swap the router URL, mirror it into location, then fire the afterNavigate callback
  // registered by bind_comparison_url
  let on_navigation: Parameters<typeof afterNavigate>[0]
  const navigate = (url: string, type: `enter` | `link`) => {
    page.url = new URL(url, location.origin)
    history.replaceState(null, ``, page.url.href)
    on_navigation({ type } as never)
    flushSync()
  }
  let cleanup = () => {}

  beforeEach(() => {
    vi.mocked(afterNavigate).mockImplementationOnce((callback) => {
      on_navigation = callback
    })
    cleanup = $effect.root(() => bind_comparison_url())
  })
  afterEach(() => cleanup()) // tear down the effect root even when an assertion failed

  it.each([
    [`/`, true],
    [`/models`, true],
    [`/models/${key_a}`, true],
    [`/models/unknown-model`, false],
    [`/models/${key_a}/extra`, false],
    [`/benchmarks/discovery`, true],
    [`/benchmarks/geo-opt`, true],
    [`/benchmarks/phonons`, true],
    [`/benchmarks/md`, true],
    [`/benchmarks/diatomics`, true],
    [`/data`, false],
    [`/data/sets`, false],
    [`/data/wbm`, false],
    [`/data/tmi`, false],
    [`/api`, false],
    [`/contribute`, false],
    [`/changelog`, false],
    [`/benchmarks`, false],
    [`/benchmarks/diatomics/tmi`, false],
    [`/benchmarks/discovery/tmi`, false],
  ])(`scopes shared comparison links on %s (enabled=%s)`, (route, enabled) => {
    navigate(`${route}?compare=${key_a},${key_b}&keep=1#section`, `enter`)

    expect(comparison.open).toBe(enabled)
    expect([...comparison.keys]).toEqual(enabled ? [key_a, key_b] : [])
    expect(location.search).toBe(
      enabled ? `?compare=${key_a},${key_b}&keep=1` : `?keep=1`,
    )
    expect(location.hash).toBe(`#section`)
  })

  it(`reads shared links, keeps the selection across in-app navigation and re-applies it`, () => {
    navigate(`/?compare=${key_a},no-such-model,${key_b}`, `enter`)
    expect([...comparison.keys]).toEqual([key_a, key_b])
    expect(comparison.open).toBe(true) // shared link with >1 model opens the dialog
    expect(location.search).toBe(`?compare=${key_a},${key_b}`) // stale key dropped

    navigate(`/data?compare=${key_c}&keep=1`, `link`)
    expect([...comparison.keys]).toEqual([key_a, key_b])
    expect(comparison.open).toBe(false)
    expect(location.search).toBe(`?keep=1`)

    navigate(`/benchmarks/md`, `link`) // absent param: selection follows the user, dialog closes
    expect([...comparison.keys]).toEqual([key_a, key_b])
    expect(comparison.open).toBe(false)
    expect(location.search).toBe(`?compare=${key_a},${key_b}`)

    // same-path navigation that only drops the query (nav link of the current page) must
    // re-apply the param too
    navigate(`/benchmarks/md`, `link`)
    expect(location.search).toBe(`?compare=${key_a},${key_b}`)

    comparison.toggle(key_b)
    flushSync()
    expect(location.search).toBe(`?compare=${key_a}`)
    comparison.keys.clear()
    flushSync()
    expect(location.search).toBe(``)

    navigate(`/?compare=${key_c}`, `enter`) // a single model is nothing to compare yet
    expect([...comparison.keys]).toEqual([key_c])
    expect(comparison.open).toBe(false)

    comparison.set([key_a, key_b])
    navigate(`/`, `enter`) // cached selection alone is not a shared link
    expect(comparison.open).toBe(false)
  })
})

describe(`ModelComparison dialog`, () => {
  with_reactive_page_url()
  // fire only callbacks registered since from_idx (earlier tests' components are gone)
  let callback_start_idx = 0
  const navigate = async (
    path: string,
    type: AfterNavigate[`type`] = `link`,
    from_idx = callback_start_idx,
  ) => {
    page.url = new URL(path, location.origin)
    history.replaceState(null, ``, path)
    for (const [callback] of vi.mocked(afterNavigate).mock.calls.slice(from_idx)) {
      callback({ type, from: null, to: { url: page.url } } as AfterNavigate)
    }
    await tick()
  }
  beforeEach(() => {
    callback_start_idx = vi.mocked(afterNavigate).mock.calls.length
  })

  it(`opens with a single selected model without looping and lists it`, async () => {
    comparison.toggle(key_a)
    mount(ModelComparison, { target: document.body })
    await tick()
    // a selection alone renders nothing outside the dialog (no floating tray)
    expect(document.querySelector(`dialog[open]`)).toBeNull()
    expect(document.body.textContent).not.toContain(MODELS[0].model_name)
    comparison.open = true
    // a few flushes: a reactive write-back loop would throw effect_update_depth_exceeded
    // here instead of settling
    for (let idx = 0; idx < 5; idx++) await tick()

    expect(comparison.open).toBe(true)
    expect([...comparison.keys]).toEqual([key_a])
    const header_cells = [...document.querySelectorAll(`thead th a`)]
    expect(header_cells.map((cell) => cell.textContent)).toEqual([MODELS[0].model_name])
    // with nothing to compare yet, the model picker is the next step and takes focus
    expect(document.activeElement).toBe(
      document.querySelector(`dialog input[role="combobox"]`),
    )
  })

  it(`the model page opens comparison without deselecting the current model`, async () => {
    mount(ModelPage, {
      target: document.body,
      props: { data: { model_key: key_a, md_per_system: null } },
    })
    const button = doc_query<HTMLButtonElement>(`.links > button`)
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
    comparison.toggle(key_c)
    await tick()
    expect(button.textContent?.trim()).toBe(`Comparing with 2 others`)
    expect(button.classList.contains(`selected`)).toBe(true)
  })

  it(`removes models via the header buttons and adds via the picker`, async () => {
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
      doc_query<HTMLButtonElement>(`thead th button[aria-label^="Remove"]`).click()
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
    const sort_select = doc_query<HTMLSelectElement>(`.option-sort select`)
    sort_select.value = `added`
    sort_select.dispatchEvent(new Event(`change`, { bubbles: true }))
    await tick()
    const newest = unselected.toSorted((md_1, md_2) =>
      (md_2.dates.benchmark_added ?? ``).localeCompare(md_1.dates.benchmark_added ?? ``),
    )[0]
    expect(option_labels()[0]).toContain(`${newest.model_name} `)
    expect(option_labels()[0]).toContain(`· ${newest.dates.benchmark_added}`)

    // picking an option adds that model (MultiSelect's bind:value writes back to the
    // store), and the dialog's own close button closes it
    doc_query<HTMLLIElement>(`dialog ul.options li`).click()
    await tick()
    expect([...comparison.keys]).toEqual([key_c, newest.model_key])
    expect(document.querySelectorAll(`thead th a`)).toHaveLength(2)
    doc_query<HTMLButtonElement>(`button[aria-label="Close model comparison"]`).click()
    await tick()
    expect(comparison.open).toBe(false)
  })

  it(`plots cost vs accuracy via DynamicScatter, y defaulting to the current task page's metric`, async () => {
    page.url = new URL(`http://localhost/benchmarks/phonons`)
    const with_kappa = MODELS.filter((model) => model.metrics?.phonons?.kappa_103?.κ_SRME)
    const compared = with_kappa.slice(0, 2).map((model) => model.model_key)
    comparison.set(compared)
    comparison.open = true
    mount(ModelComparison, { target: document.body })
    await tick()
    const plot = () => get_scatter_plot_props(plot_mocks.ScatterPlot) as PlotProps
    expect(plot().x_axis).toMatchObject({
      label: `Params`,
      scale_type: `log`,
    })
    expect(plot().y_axis).toMatchObject({
      label: `κ<sub>SRME</sub>`,
    })
    // one series per model with the dimmed field first, highlighted (ringed, labeled)
    // compared models last, then the dashed Pareto frontier
    const series = plot().series
    const frontier = series.at(-1)
    expect(frontier?.markers).toBe(`line`)
    // only models with finite x, y, color and size values get a point
    const model_series = series.slice(0, -1)
    expect(model_series.length).toBeGreaterThan(compared.length)
    expect(model_series.length).toBeLessThanOrEqual(ACTIVE_MODELS.length)
    const highlighted = model_series.filter((trace) => trace.point_style?.stroke)
    expect(highlighted.map((trace) => trace.id)).toEqual(compared)
    expect(model_series.slice(-2)).toEqual(highlighted)
    expect(highlighted.every((trace) => trace.point_label?.text)).toBe(true)
    const dimmed = model_series.slice(0, -2)
    expect(dimmed.every((trace) => trace.point_style?.fill_opacity === 0.3)).toBe(true)
    expect(dimmed.every((trace) => trace.point_label === undefined)).toBe(true)

    // Axis selections update the derived series; clicking a point toggles the model.
    plot().on_axis_change(`x`, `n_training_materials`)
    plot().on_axis_change(`y`, `F1`)
    await tick()
    expect(plot().x_axis.label).toBe(`N<sub>materials</sub>`)
    expect(plot().y_axis.label).toBe(`F1`)
    expect(query_param(`compare_plot_x`)).toBe(`n_training_materials`)
    comparison.open = false
    await tick()
    comparison.open = true
    await tick()
    expect(plot().x_axis.label).toBe(`N<sub>materials</sub>`)

    // A retained dialog follows the new task's default, while explicit shared axes win.
    await navigate(`/benchmarks/md?compare_plot_x=n_training_materials`)
    expect(plot().y_axis.label).toBe(`CMDS`)
    expect(query_param(`compare_plot_y`)).toBeNull()
    expect(query_param(`compare_plot_x`)).toBe(`n_training_materials`)

    // Emptying the selection destroys the scatter; rebuilding it must not turn the
    // current x selection into a new default and silently erase it from the URL.
    comparison.keys.clear()
    await tick()
    const callback_idx = vi.mocked(afterNavigate).mock.calls.length
    comparison.set(compared)
    await tick()
    await navigate(`${page.url.pathname}${page.url.search}`, `enter`, callback_idx)
    expect(plot().x_axis.label).toBe(`N<sub>materials</sub>`)
    expect(query_param(`compare_plot_x`)).toBe(`n_training_materials`)
    plot().point_events.onclick({ point: { metadata: { model_key: compared[0] } } })
    expect([...comparison.keys]).toEqual([compared[1]])
  })

  it(`restores a late-mounted plot's shared linear axes and preserves them across score changes`, async () => {
    const cleanup = $effect.root(() => bind_score_weights())
    try {
      await navigate(
        `/?cps_weights=1,0,0&compare_plot_y=model_params&compare_plot_x_scale=linear&compare_plot_y_scale=linear`,
      )
      comparison.set([key_a, key_b])
      comparison.open = true
      // Kit registers late components for future navigations without replaying the
      // initial event. Suppress tracked_mount's eager replay to match that lifecycle.
      vi.mocked(afterNavigate).mockImplementationOnce(() => {})
      mount(ModelComparison, { target: document.body })
      await tick()
      const plot = () => get_scatter_plot_props(plot_mocks.ScatterPlot) as PlotProps
      expect(score_weight_records().CPS).toEqual({ F1: 1, κ_SRME: 0, RMSD: 0 })
      for (const axis of [`x`, `y`] as const) {
        expect(plot()[`${axis}_axis`]).toMatchObject({
          label: `Params`,
          scale_type: `linear`,
        })
        expect(query_param(`compare_plot_${axis}_scale`)).toBe(`linear`)
      }
      await navigate(
        `/?cps_weights=0,1,0&compare_plot_y=model_params&compare_plot_x_scale=linear&compare_plot_y_scale=linear`,
      )
      expect(score_weight_records().CPS).toEqual({ F1: 0, κ_SRME: 1, RMSD: 0 })
      expect([plot().x_axis.scale_type, plot().y_axis.scale_type]).toEqual([
        `linear`,
        `linear`,
      ])
      expect(query_param(`cps_weights`)).toBe(`0,1,0`)
    } finally {
      await navigate(`/`, `popstate`)
      cleanup()
    }
  })
})
