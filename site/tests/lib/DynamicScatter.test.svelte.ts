import { HYPERPARAMS, METADATA_COLS } from '$lib/labels'
import DynamicScatter from '$lib/plot/DynamicScatter.svelte'
import type { ModelData } from '$lib/types'
import { tick } from 'svelte'
import { expect, it } from 'vitest'
import { mount } from '../index'

const make_models = (min_value: number, max_value: number): ModelData[] =>
  [min_value, max_value].map(
    (value) =>
      ({
        model_key: `model-${value}`,
        model_name: `Model ${value}`,
        date_added: `2025-01-01`,
        model_params: value,
        n_training_materials: value,
        n_training_structures: value,
      }) as ModelData,
  )

it.each([
  {
    scenario: `99-fold positive range`,
    min_value: 1,
    max_value: 99,
    x_key: HYPERPARAMS.model_params.key,
    expected: [false, false, false, false],
  },
  {
    scenario: `100-fold positive range`,
    min_value: 1,
    max_value: 100,
    x_key: HYPERPARAMS.model_params.key,
    expected: [true, true, true, true],
  },
  {
    scenario: `non-positive minimum`,
    min_value: 0,
    max_value: 100,
    x_key: HYPERPARAMS.model_params.key,
    expected: [false, false, false, false],
  },
  {
    scenario: `date x-axis`,
    min_value: 1,
    max_value: 100,
    x_key: METADATA_COLS.date_added.key,
    expected: [false, true, true, true],
  },
])(`sets log toggles for $scenario`, ({ min_value, max_value, x_key, expected }) => {
  mount(DynamicScatter, {
    target: document.body,
    props: {
      models: make_models(min_value, max_value),
      x_key,
      y_key: METADATA_COLS.n_training_materials.key,
      color_key: METADATA_COLS.n_training_structures.key,
      show_model_labels: false,
    },
  })

  const toggles = [...document.querySelectorAll<HTMLInputElement>(`.log-controls input`)]
  expect(toggles.map(({ checked, disabled }) => ({ checked, disabled }))).toStrictEqual(
    expected.map((enabled) => ({ checked: enabled, disabled: !enabled })),
  )
})

it(`re-evaluates manual log choices after an axis change`, async () => {
  let x_key = $state(HYPERPARAMS.model_params.key)
  mount(DynamicScatter, {
    target: document.body,
    props: {
      models: make_models(1, 100),
      get x_key() {
        return x_key
      },
      set x_key(value) {
        x_key = value
      },
      y_key: METADATA_COLS.n_training_materials.key,
      color_key: METADATA_COLS.n_training_structures.key,
      show_model_labels: false,
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
