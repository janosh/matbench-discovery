import { goto } from '$app/navigation'
import { HYPERPARAMS, METADATA_COLS } from '$lib/labels'
import DynamicScatter from '$lib/plot/DynamicScatter.svelte'
import { interpolateViridis } from 'd3-scale-chromatic'
import { tick } from 'svelte'
import { beforeEach, expect, it, vi } from 'vitest'
import {
  choose_scatter_property,
  doc_query,
  mount,
  mount_with_url,
  navigate,
  query_param,
} from '../index'

const make_models = (...values: number[]) =>
  values.map((value) => ({
    model_key: `model-${value}`,
    model_name: `Model ${value}`,
    dates: { benchmark_added: `2025-01-01` },
    model_params: value,
    n_training_materials: value,
    n_training_structures: value,
  }))
const scatter_props = {
  y_key: METADATA_COLS.n_training_materials.key,
  color_key: METADATA_COLS.n_training_structures.key,
  show_model_labels: false,
}
const marker_color = (marker: Element) =>
  marker
    .closest<SVGElement>(`[style*="--point-fill-color"]`)
    ?.style.getPropertyValue(`--point-fill-color`)

beforeEach(() => {
  vi.spyOn(HTMLElement.prototype, `clientWidth`, `get`).mockReturnValue(800)
  vi.spyOn(HTMLElement.prototype, `clientHeight`, `get`).mockReturnValue(600)
})

it.each([
  [1, 99, HYPERPARAMS.model_params.key, []],
  [1, 100, HYPERPARAMS.model_params.key, [`X`, `Y`, `Color`, `Size`]],
  [0, 100, HYPERPARAMS.model_params.key, []],
  [1, 100, METADATA_COLS.benchmark_added.key, [`Y`, `Color`, `Size`]],
] as const)(
  `sets log toggles for range [%d, %d] with x=%s`,
  async (min_value, max_value, x_key, expected_labels) => {
    mount(DynamicScatter, {
      target: document.body,
      props: {
        models: make_models(min_value, max_value),
        x_key,
        ...scatter_props,
      },
    })
    await tick()

    const toggles = [
      ...document.querySelectorAll<HTMLInputElement>(`.log-controls input`),
    ]
    expect(
      toggles.map((toggle) => toggle.parentElement?.textContent?.trim()),
    ).toStrictEqual(expected_labels)
    expect(toggles.every((toggle) => toggle.checked && !toggle.disabled)).toBe(true)
    expect(document.querySelector(`.log-controls`) !== null).toBe(
      expected_labels.length > 0,
    )
    const plot_area = doc_query<SVGRectElement>(`.scatter clipPath rect`)
    const label_group = doc_query(`.scatter .y-label`).closest(`g`)
    const plot_center_y =
      Number(plot_area.getAttribute(`y`)) + Number(plot_area.getAttribute(`height`)) / 2
    const rotation = label_group?.getAttribute(`transform`)
    expect(rotation).toMatch(/^rotate\(-90, /)
    expect(Number(rotation?.split(`, `)[2]?.replace(`)`, ``))).toBe(plot_center_y)
  },
)

it(`re-evaluates manual log choices after an axis change`, async () => {
  const axis = $state({ key: HYPERPARAMS.model_params.key })
  mount(DynamicScatter, {
    target: document.body,
    props: {
      models: make_models(1, 100),
      ...scatter_props,
      get x_key() {
        return axis.key
      },
      set x_key(value) {
        axis.key = value
      },
    },
  })
  await tick()

  const x_ticks = () =>
    [...document.querySelectorAll(`.x-axis .tick text`)].map((label) => label.textContent)
  const log_ticks = x_ticks()
  const x_toggle = document.querySelector<HTMLInputElement>(`.log-controls input`)
  expect(x_toggle?.checked).toBe(true)
  x_toggle?.click()
  await tick()
  expect(x_toggle?.checked).toBe(false)
  expect(x_ticks()).not.toEqual(log_ticks)
  expect(query_param(`x_scale`)).toBe(`linear`)
  const shared_url = location.href
  const linear_ticks = x_ticks()
  await navigate(`http://localhost/`, `popstate`)
  expect(x_toggle?.checked).toBe(true)
  await navigate(shared_url, `popstate`)
  expect(x_toggle?.checked).toBe(false)
  expect(x_ticks()).toEqual(linear_ticks)

  await choose_scatter_property(`X axis`, `Training Materials`)
  expect(axis.key).toBe(METADATA_COLS.n_training_materials.key)
  expect(x_toggle?.checked).toBe(true)
  expect(x_ticks()).toEqual(log_ticks)
  expect(query_param(`x_scale`)).toBeNull()
  await choose_scatter_property(`X axis`, `Params`)
  expect(x_toggle?.checked).toBe(true)
  expect(query_param(`x_scale`)).toBeNull()
})

it.each([
  { counts: [0, 1, 10, 40], scale: `Arcsinh`, fraction: Math.asinh(4) / Math.asinh(160) },
  { counts: [1, 2, 10, 40], scale: `Logarithmic`, fraction: Math.log(2) / Math.log(40) },
  { counts: [0, 0, 0, 0], scale: `Arcsinh`, fraction: 0.5 },
])(
  `honors log preferences with $scale color scaling for $counts`,
  async ({ counts, scale, fraction }) => {
    mount(DynamicScatter, {
      target: document.body,
      props: {
        models: counts.map((count, idx) => ({
          ...make_models(idx + 1)[0],
          n_training_structures: count,
        })),
        options: [
          { ...HYPERPARAMS.model_params },
          { ...METADATA_COLS.n_training_materials },
          { ...METADATA_COLS.n_training_structures, scale_type: `log` },
        ],
        x_key: HYPERPARAMS.model_params.key,
        ...scatter_props,
      },
    })
    await tick()
    const toggle = doc_query<HTMLInputElement>(`[aria-label="${scale} scales"] input`)
    expect(toggle.parentElement?.textContent?.trim()).toBe(`Color`)
    expect(toggle.checked).toBe(true)
    const colors = () => [...document.querySelectorAll(`path.marker`)].map(marker_color)
    expect(colors()).toHaveLength(counts.length)
    expect(colors()[1]).toBe(interpolateViridis(fraction))
    toggle.click()
    await tick()
    expect(toggle.checked).toBe(false)
    const [min, value, , max] = counts
    expect(colors()[1]).toBe(
      interpolateViridis(max === min ? 0.5 : (value - min) / (max - min)),
    )
    expect(colors()).toHaveLength(counts.length)
  },
)

it.each([
  [`x`, `X axis`],
  [`y`, `Y axis`],
  [`color`, `Color`],
  [`size`, `Marker size`],
] as const)(
  `selects the %s property and updates its data binding and rendered plot`,
  async (dim, label) => {
    const models = [1, 4, 10].map((value, idx) => ({
      ...make_models(value)[0],
      n_training_materials: idx + 1,
    }))
    const selection = $state({ key: `model_params` })
    const prop_key = `${dim}_key` as const
    mount(DynamicScatter, {
      target: document.body,
      props: {
        models,
        show_model_labels: false,
        x_key: `model_params`,
        y_key: `model_params`,
        color_key: `n_training_structures`,
        size_key: `model_params`,
        get [prop_key]() {
          return selection.key
        },
        set [prop_key](value: string) {
          selection.key = value
        },
      },
    })
    await tick()
    expect(
      document.querySelectorAll(`.property-picker input[role="combobox"]`),
    ).toHaveLength(1)
    expect(document.querySelectorAll(`.interactive-axis-label button`)).toHaveLength(2)
    const rendered_values = () =>
      dim === `color` || dim === `size`
        ? [...document.querySelectorAll(`path.marker`)].map((marker) =>
            dim === `color` ? marker_color(marker) : marker.getAttribute(`d`),
          )
        : [...document.querySelectorAll(`.${dim}-axis .tick text`)].map(
            (tick_label) => tick_label.textContent,
          )
    const previous_values = rendered_values()
    expect(previous_values.length).toBeGreaterThan(0)
    await choose_scatter_property(label, `Training Materials`)
    expect(selection.key).toBe(`n_training_materials`)
    expect(rendered_values()).not.toEqual(previous_values)
    expect(query_param(dim)).toBe(`n_training_materials`)
    const shared_url = location.href
    const selected_values = rendered_values()
    await navigate(`http://localhost/`, `popstate`)
    expect(selection.key).toBe(`model_params`)
    expect(rendered_values()).toEqual(previous_values)
    await navigate(shared_url, `popstate`)
    expect(selection.key).toBe(`n_training_materials`)
    expect(rendered_values()).toEqual(selected_values)
  },
)

it(`keeps duplicate-label series distinct and collapses their legend`, async () => {
  const models = make_models(10, 20).map((model) => ({
    ...model,
    model_name: `Duplicate label`,
  }))
  mount(DynamicScatter, {
    target: document.body,
    props: {
      models,
      x_key: METADATA_COLS.benchmark_added.key,
      ...scatter_props,
    },
  })

  doc_query(`button.models-toggle`).click()
  await tick()
  const legend = doc_query(`.scatter > .legend:has(.legend-item)`)
  expect(
    [...legend.querySelectorAll(`.legend-item`)].map((item) => item.textContent?.trim()),
  ).toEqual(models.map(({ model_name, model_key }) => `${model_name} (${model_key})`))
  expect(document.querySelector(`button.models-toggle`)).toBeNull()

  const click = () => new MouseEvent(`click`, { bubbles: true })
  doc_query(`.legend-item`, legend).dispatchEvent(click())
  await tick()
  expect(document.querySelector(`button.models-toggle`)).toBeNull()

  // the controls row sits inside the same wrapper the attachment is on, so only the
  // legend itself may count as inside — everything else collapses it
  doc_query(`.collapsible-legend .controls-row`).dispatchEvent(click())
  await tick()
  expect(document.querySelector(`button.models-toggle`)).not.toBeNull()
})

it(`isolates comparison selections and linear scales from the surrounding page URL`, async () => {
  const page_params = {
    x: `n_training_materials`,
    y: `n_training_structures`,
    color: `model_params`,
    size: `n_training_materials`,
  }
  const comparison_params = {
    compare_plot_x: `n_training_structures`,
    compare_plot_y: `model_params`,
    compare_plot_color: `n_training_materials`,
    compare_plot_size: `n_training_structures`,
    ...Object.fromEntries(
      [`x`, `y`, `color`, `size`].map((dim) => [`compare_plot_${dim}_scale`, `linear`]),
    ),
  }
  await mount_with_url(
    DynamicScatter,
    `http://localhost/?${new URLSearchParams({ ...page_params, ...comparison_params })}`,
    {
      props: {
        models: make_models(1, 100),
        x_key: `model_params`,
        ...scatter_props,
        url_prefix: `compare_plot`,
      },
    },
  )
  expect(doc_query(`.x-label`).textContent).toContain(`Training Structures`)
  expect(doc_query(`.y-label`).textContent).toContain(`Params`)
  expect(doc_query(`.colorbar .property-select`).textContent).toContain(
    `Training Materials`,
  )
  expect(doc_query(`.property-picker .selected-label`).textContent).toContain(
    `Training Structures`,
  )
  const toggles = [...document.querySelectorAll<HTMLInputElement>(`.log-controls input`)]
  expect(toggles).toHaveLength(4)
  expect(toggles.every(({ checked }) => !checked)).toBe(true)
  for (const [key, value] of Object.entries(comparison_params))
    expect(query_param(key)).toBe(value)

  await choose_scatter_property(`X axis`, `Training Materials`)
  expect(query_param(`compare_plot_x`)).toBe(`n_training_materials`)
  expect(query_param(`compare_plot_x_scale`)).toBeNull()
  for (const [key, value] of Object.entries(page_params))
    expect(query_param(key)).toBe(value)
})

it(`rejects categorical axes and sizes while allowing categorical colors`, async () => {
  await mount_with_url(
    DynamicScatter,
    `http://localhost/?x=category&y=missing&size=category&color=category&x_scale=invalid`,
    {
      props: {
        models: make_models(1, 100).map((model) => ({ ...model, category: `public` })),
        x_key: `model_params`,
        y_key: `n_training_materials`,
        color_key: `category`,
        options: [
          HYPERPARAMS.model_params,
          METADATA_COLS.n_training_materials,
          { key: `category`, label: `Category`, categories: { public: `red` } },
        ],
      },
    },
  )
  expect(doc_query(`.x-label`).textContent).toContain(`Params`)
  expect(doc_query(`.y-label`).textContent).toContain(`Training Materials`)
  expect(doc_query(`.colorbar .property-select`).textContent).toContain(`Category`)
  expect(doc_query(`.property-picker .selected-label`).textContent).toContain(`Params`)
  expect(location.search).toBe(``)
})

it(`uses the active task's defaults when a comparison plot survives navigation`, async () => {
  const defaults = $state({ x: `model_params`, y: `n_training_materials` })
  await mount_with_url(
    DynamicScatter,
    `http://localhost/benchmarks/md?compare_plot_x=n_training_structures`,
    {
      props: {
        models: make_models(1, 100),
        x_key: `n_training_structures`,
        ...scatter_props,
        url_prefix: `compare_plot`,
        get url_defaults() {
          return defaults
        },
      },
    },
  )
  expect(query_param(`compare_plot_x`)).toBe(`n_training_structures`)
  await choose_scatter_property(`Y axis`, `Params`)
  expect(query_param(`compare_plot_y`)).toBe(`model_params`)
  defaults.y = `n_training_structures`
  await navigate(`http://localhost/benchmarks/diatomics`, `popstate`)
  expect(doc_query(`.x-label`).textContent).toContain(`Params`)
  expect(doc_query(`.y-label`).textContent).toContain(`Training Structures`)
  expect(query_param(`compare_plot_y`)).toBeNull()
  await navigate(
    `http://localhost/benchmarks/diatomics?compare_plot_y=model_params`,
    `popstate`,
  )
  expect(doc_query(`.y-label`).textContent).toContain(`Params`)
  expect(query_param(`compare_plot_y`)).toBe(`model_params`)
})

it(`renders category colors and dataset links without model metadata`, async () => {
  const datasets = [
    { slug: `public-set`, name: `Public set`, count: 100, access: `public` },
    { slug: `partial-set`, name: `Partial set`, count: 200, access: `partial` },
    { slug: `unknown-access`, name: `Unknown access`, count: 300, access: `unknown` },
    { slug: `missing-count`, name: `Missing count`, count: null, access: `public` },
  ]
  const categories = { public: `#25836d`, partial: `#b16c00` }
  const props = {
    models: datasets,
    get_identity: ({ slug, name }: (typeof datasets)[number]) => ({
      key: slug,
      name,
      href: `/data/${slug}`,
    }),
    options: [
      { key: `count`, label: `Count`, description: `Structure count` },
      { key: `unreported`, label: `Unreported`, description: `No values available` },
      {
        key: `access`,
        label: `Access`,
        description: `Corpus availability`,
        categories,
      },
    ],
    x_key: `count`,
    y_key: `count`,
    color_key: `access`,
    size_key: `count`,
    legend: null,
    hovered: true,
  }
  mount<typeof props, Record<string, unknown>>(DynamicScatter, {
    target: document.body,
    props,
  })
  await tick()
  const markers = [...document.querySelectorAll<SVGPathElement>(`path.marker`)]
  expect(markers).toHaveLength(2)
  expect(doc_query(`.property-picker .selected-label small`).textContent).toMatch(
    /3\s+models/,
  )
  const category_fills = markers.map(marker_color)
  expect(doc_query(`.colorbar .property-select`).textContent).toContain(`Access`)
  expect(doc_query(`.colorbar-wrapper`).getAttribute(`role`)).toBe(`group`)
  expect(document.querySelector(`.colorbar .bar`)).toBeNull()
  for (const [idx, color] of Object.values(categories).entries()) {
    expect(category_fills[idx]).toBe(color)
  }
  markers[1].dispatchEvent(new MouseEvent(`click`, { bubbles: true }))
  expect(goto).toHaveBeenCalledWith(`/data/partial-set`)
  await choose_scatter_property(`Color`, `Count`)
  expect(document.querySelector(`.category-legend`)).toBeNull()
  expect(document.querySelector(`.colorbar .bar`)).not.toBeNull()
  expect(document.querySelectorAll(`path.marker`)).toHaveLength(3)
  expect(markers.map(marker_color)).not.toEqual(category_fills)
  doc_query<SVGPathElement>(`path.marker`).dispatchEvent(
    new MouseEvent(`click`, { bubbles: true }),
  )
  await tick()
  await choose_scatter_property(`Color`, `Access`)
  expect(document.querySelectorAll(`path.marker`)).toHaveLength(2)
  await choose_scatter_property(`Color`, `Unreported`)
  expect(document.querySelectorAll(`path.marker`)).toHaveLength(0)
  expect(document.querySelector(`.colorbar .bar`)).toBeNull()
  expect(doc_query(`.colorbar .property-select`).textContent).toContain(`Unreported`)
  await choose_scatter_property(`Color`, `Access`)
  expect(document.querySelectorAll(`path.marker`)).toHaveLength(2)
  expect([...document.querySelectorAll(`path.marker`)].map(marker_color)).toEqual(
    category_fills,
  )
})

it(`dims and unlabels models outside highlight_keys, drawing highlighted ones last`, async () => {
  const models = make_models(1, 100)
  mount(DynamicScatter, {
    target: document.body,
    props: {
      models,
      x_key: HYPERPARAMS.model_params.key,
      highlight_keys: new Set([models[0].model_key]),
      ...scatter_props,
      show_model_labels: true,
      bleed: false,
      legend: null,
    },
  })
  await tick()

  const markers = [...document.querySelectorAll<SVGPathElement>(`path.marker`)]
  const style_of = (marker: SVGPathElement) => ({
    fill_opacity: marker.getAttribute(`fill-opacity`),
    stroke: marker.getAttribute(`stroke`),
  })
  // the dimmed model paints first (default stroke), the ringed highlighted one on top
  expect(markers.map(style_of)).toEqual([
    { fill_opacity: `0.3`, stroke: `#000` },
    { fill_opacity: `1`, stroke: `currentColor` },
  ])
  const labels = [...document.querySelectorAll(`.scatter text`)]
    .map((text) => text.textContent?.trim())
    .filter((text) => text?.startsWith(`Model `))
  expect(labels).toEqual([models[0].model_name])
  // no full-bleed wrapper and no collapsed-legend toggle inside dialogs
  expect(document.querySelector(`.bleed-1400`)).toBeNull()
  expect(document.querySelector(`button.models-toggle`)).toBeNull()
})
