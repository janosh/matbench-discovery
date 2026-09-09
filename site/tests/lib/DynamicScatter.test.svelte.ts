import { HYPERPARAMS, METADATA_COLS } from '$lib/labels'
import DynamicScatter from '$lib/plot/DynamicScatter.svelte'
import type { ModelData } from '$lib/types'
import { tick } from 'svelte'
import { beforeEach, expect, it, vi } from 'vitest'
import { choose_scatter_property, doc_query, mount } from '../index'

const make_models = (min_value: number, max_value: number): ModelData[] =>
  [min_value, max_value].map(
    (value) =>
      ({
        model_key: `model-${value}`,
        model_name: `Model ${value}`,
        dates: { benchmark_added: `2025-01-01` },
        model_params: value,
        n_training_materials: value,
        n_training_structures: value,
      }) as ModelData,
  )
const scatter_props = {
  y_key: METADATA_COLS.n_training_materials.key,
  color_key: METADATA_COLS.n_training_structures.key,
  show_model_labels: false,
}

beforeEach(() => {
  vi.spyOn(HTMLElement.prototype, `clientWidth`, `get`).mockReturnValue(800)
  vi.spyOn(HTMLElement.prototype, `clientHeight`, `get`).mockReturnValue(600)
})

it.each([
  {
    scenario: `99-fold positive range`,
    min_value: 1,
    max_value: 99,
    x_key: HYPERPARAMS.model_params.key,
    expected_labels: [],
  },
  {
    scenario: `100-fold positive range`,
    min_value: 1,
    max_value: 100,
    x_key: HYPERPARAMS.model_params.key,
    expected_labels: [`X`, `Y`, `Color`, `Size`],
  },
  {
    scenario: `non-positive minimum`,
    min_value: 0,
    max_value: 100,
    x_key: HYPERPARAMS.model_params.key,
    expected_labels: [],
  },
  {
    scenario: `date x-axis`,
    min_value: 1,
    max_value: 100,
    x_key: METADATA_COLS.benchmark_added.key,
    expected_labels: [`Y`, `Color`, `Size`],
  },
])(
  `sets log toggles for $scenario`,
  async ({ min_value, max_value, x_key, expected_labels }) => {
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

  await choose_scatter_property(`X axis`, `Training Materials`)
  expect(axis.key).toBe(METADATA_COLS.n_training_materials.key)
  expect(x_toggle?.checked).toBe(true)
  expect(x_ticks()).toEqual(log_ticks)
})

it.each([
  [`x`, `X axis`],
  [`y`, `Y axis`],
  [`color`, `Color`],
  [`size`, `Marker size`],
] as const)(
  `searches the %s picker and updates its data binding and rendered plot`,
  async (dim, label) => {
    const models = [1, 4, 10].map((value, idx) => ({
      ...make_models(value, value)[0],
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
    ).toHaveLength(4)
    const rendered_values = () =>
      dim === `color` || dim === `size`
        ? [...document.querySelectorAll(`path.marker`)].map((marker) =>
            marker.getAttribute(dim === `color` ? `fill` : `d`),
          )
        : [...document.querySelectorAll(`.${dim}-axis .tick text`)].map(
            (tick_label) => tick_label.textContent,
          )
    const previous_values = rendered_values()
    expect(previous_values.length).toBeGreaterThan(0)
    await choose_scatter_property(label, `Training Materials`)
    expect(selection.key).toBe(`n_training_materials`)
    expect(rendered_values()).not.toEqual(previous_values)
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
