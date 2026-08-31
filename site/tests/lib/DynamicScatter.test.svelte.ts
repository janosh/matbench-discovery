import { HYPERPARAMS, METADATA_COLS, scatter_options } from '$lib/labels'
import DynamicScatter from '$lib/plot/DynamicScatter.svelte'
import type { ModelData } from '$lib/types'
import { tick } from 'svelte'
import { expect, it, vi } from 'vitest'
import { doc_query, mount } from '../index'

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
  ({ min_value, max_value, x_key, expected_labels }) => {
    mount(DynamicScatter, {
      target: document.body,
      props: {
        models: make_models(min_value, max_value),
        x_key,
        ...scatter_props,
      },
    })

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
  },
)

it(`re-evaluates manual log choices after an axis change`, async () => {
  let x_key = $state(HYPERPARAMS.model_params.key)
  mount(DynamicScatter, {
    target: document.body,
    props: {
      models: make_models(1, 100),
      ...scatter_props,
      get x_key() {
        return x_key
      },
      set x_key(value) {
        x_key = value
      },
    },
  })

  const x_toggle = document.querySelector<HTMLInputElement>(`.log-controls input`)
  expect(x_toggle?.checked).toBe(true)
  x_toggle?.click()
  expect(x_toggle?.checked).toBe(false)

  x_key = METADATA_COLS.n_training_materials.key
  await tick()
  expect(x_toggle?.checked).toBe(true)
})

it.each([
  { scenario: `default threshold`, select_filter_threshold: undefined, has_filter: true },
  {
    scenario: `threshold equal to the option count`,
    select_filter_threshold: scatter_options.length,
    has_filter: false,
  },
])(
  `adds a dropdown filter for option lists above the $scenario`,
  async ({ select_filter_threshold, has_filter }) => {
    vi.spyOn(HTMLElement.prototype, `clientWidth`, `get`).mockReturnValue(800)
    vi.spyOn(HTMLElement.prototype, `clientHeight`, `get`).mockReturnValue(600)
    mount(DynamicScatter, {
      target: document.body,
      props: {
        models: make_models(1, 100),
        x_key: HYPERPARAMS.model_params.key,
        select_filter_threshold,
        ...scatter_props,
      },
    })
    await tick()

    doc_query<HTMLButtonElement>(`.axis-trigger`).click()
    await vi.waitFor(() =>
      expect(document.querySelector(`.portal-select-filter input`) !== null).toBe(
        has_filter,
      ),
    )
    const filter_input = document.querySelector<HTMLInputElement>(
      `.portal-select-filter input`,
    )

    if (filter_input) {
      filter_input.value = `no option has this label`
      filter_input.dispatchEvent(new InputEvent(`input`, { bubbles: true }))
      expect(
        document.querySelectorAll(`.portal-select-dropdown [role="option"]`),
      ).toHaveLength(0)
    }
  },
)

it(`keeps duplicate-label series distinct and collapses their legend`, async () => {
  vi.spyOn(HTMLElement.prototype, `clientWidth`, `get`).mockReturnValue(800)
  vi.spyOn(HTMLElement.prototype, `clientHeight`, `get`).mockReturnValue(600)
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
  vi.spyOn(HTMLElement.prototype, `clientWidth`, `get`).mockReturnValue(800)
  vi.spyOn(HTMLElement.prototype, `clientHeight`, `get`).mockReturnValue(600)
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
