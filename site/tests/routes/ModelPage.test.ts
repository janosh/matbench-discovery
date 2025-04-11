import { MODELS, type ModelData } from '$lib'
import { get_metric_value } from '$lib/metrics'
import ModelPage from '$routes/models/[slug]/+page.svelte'
import { mount } from 'svelte'
import { describe, expect, it } from 'vitest'

const model_keys = MODELS.map((m) => m.model_key)
const model_key = model_keys[0]
const test_model = MODELS.find((m) => m.model_key === model_key)
if (!test_model) throw `missing test model`

describe(`Model Detail Page`, () => {
  it(`renders model details correctly`, async () => {
    // First verify test model exists
    mount(ModelPage, { target: document.body, props: { data: { model: test_model } } })

    // Check basic model info
    expect(document.querySelector(`h1`)?.textContent).toBe(test_model.model_name)
    expect(document.body.textContent).toContain(test_model.model_version)
    expect(document.body.textContent).toContain(test_model.date_added)
    expect(document.body.textContent).toContain(test_model.date_published)

    // Check meta info section
    const meta_info = document.querySelector(`.meta-info`)
    expect(meta_info?.textContent).toContain(`parameters`)
    if (test_model.n_estimators > 1) {
      expect(meta_info?.textContent).toContain(
        `Ensemble ${test_model.n_estimators} models`,
      )
    }

    // Check links section
    const links = document.querySelectorAll(`.links a`)
    const expected_link_count = [
      test_model.repo,
      test_model.paper,
      test_model.url,
      test_model.doi,
      test_model.pypi,
    ].filter(Boolean).length
    expect(links.length).toBeGreaterThanOrEqual(expected_link_count)

    // Check authors section
    const authors = document.querySelectorAll(`.authors li`)
    expect(authors).toHaveLength(test_model.authors.length)
    for (const [idx, yaml_author] of test_model.authors.entries()) {
      const author_elem = authors[idx]
      expect(author_elem.textContent).toContain(yaml_author.name)
      if (yaml_author.affiliation) {
        expect(author_elem.textContent).toContain(yaml_author.affiliation)
      }
      if (yaml_author.email) {
        expect(author_elem.querySelector(`[href^="mailto:"]`)).toBeTruthy()
      }
      if (yaml_author.github) {
        expect(author_elem.querySelector(`[href="${yaml_author.github}"]`)).toBeTruthy()
      }
      if (yaml_author.orcid) {
        expect(author_elem.querySelector(`[href="${yaml_author.orcid}"]`)).toBeTruthy()
      }
    }

    // Check trained by section if present
    if (test_model.trained_by) {
      const trainers = document.querySelectorAll(`.trained-by li`)
      expect(trainers).toHaveLength(test_model.trained_by.length)
      for (const [idx, trainer] of test_model.trained_by.entries()) {
        const trainer_el = trainers[idx]
        expect(trainer_el.textContent).toContain(trainer.name)
        if (trainer.affiliation) {
          expect(trainer_el.textContent).toContain(trainer.affiliation)
        }
      }
    }

    // Check model info section
    const model_info = document.querySelector(`.model-info`)
    expect(model_info?.textContent).toContain(test_model.model_type)
    expect(model_info?.textContent).toContain(test_model.targets)
    expect(model_info?.textContent).toContain(test_model.openness)

    // Check training set section if present
    if (test_model.training_set) {
      const training_set = document.querySelector(`.training-set`)
      expect(training_set).toBeTruthy()
      const train_set_links = training_set?.querySelectorAll(`a`)
      expect(train_set_links?.length).toBe(
        Array.isArray(test_model.training_set) ? test_model.training_set.length : 1,
      )
    }

    // Check hyperparameters section if present
    if (test_model.hyperparams) {
      const hyperparams = document.querySelector(`.hyperparams`)
      expect(hyperparams).toBeTruthy()
      for (const [key, value] of Object.entries(test_model.hyperparams)) {
        expect(hyperparams?.textContent).toContain(key)
        expect(hyperparams?.textContent).toContain(JSON.stringify(value))
      }
    }
  })
})

describe(`Models Sorting Logic`, () => {
  // Extract the common function to reduce duplication
  function sort_models(
    models: ModelData[],
    sort_by: string,
    order: `asc` | `desc`,
  ): ModelData[] {
    const sort_factor = order === `asc` ? -1 : 1

    return [...models].sort((model_1, model_2) => {
      // Special case for model_name sorting
      if (sort_by === `model_name`) {
        return order === `asc`
          ? model_1.model_name.localeCompare(model_2.model_name)
          : model_2.model_name.localeCompare(model_1.model_name)
      }

      // Get values using the helper function
      const val_1 = get_metric_value(model_1, sort_by)
      const val_2 = get_metric_value(model_2, sort_by)

      // Handle null, undefined, or NaN values by sorting last
      if (val_1 == null && val_2 == null) return 0
      if (val_1 == null || Number.isNaN(val_1)) return 1 // Always sort nulls/NaN to the end
      if (val_2 == null || Number.isNaN(val_2)) return -1 // Always sort nulls/NaN to the end

      if (typeof val_1 === `number` && typeof val_2 === `number`) {
        return sort_factor * (val_2 - val_1)
      }

      return 0
    })
  }

  // Create a set of test models that can be reused across multiple test cases
  const create_test_models = () => {
    return [
      {
        model_name: `AAA Model`,
        model_key: `aaa_model`,
        metrics: {
          discovery: {
            full_test_set: {
              F1: 0.9,
              Accuracy: 0.85,
            },
            pred_col: `is_stable`,
          },
          phonons: {
            kappa_103: { κ_SRME: 0.9 },
          },
        },
      },
      {
        model_name: `MMM Model`,
        model_key: `mmm_model`,
        metrics: {
          discovery: {
            full_test_set: {
              F1: 0.7,
              Accuracy: NaN, // Test NaN handling
            },
            pred_col: `is_stable`,
          },
          phonons: {
            kappa_103: { κ_SRME: 0.5 },
          },
        },
      },
      {
        model_name: `ZZZ Model`,
        model_key: `zzz_model`,
        metrics: {
          discovery: {
            full_test_set: {
              F1: 0.5,
              Accuracy: 0.6,
            },
            pred_col: `is_stable`,
          },
          phonons: {
            kappa_103: { κ_SRME: 0.2 },
          },
        },
      },
      {
        model_name: `Missing Data Model`,
        model_key: `missing_model`,
        metrics: {
          discovery: {
            full_test_set: { F1: 0.4 },
          },
          phonons: `not applicable` as const, // String instead of object
        },
      },
    ] as unknown as ModelData[]
  }

  // Create test models with edge cases
  const create_edge_case_models = () => {
    return [
      // Model with completely missing metrics
      {
        model_name: `No Metrics Model`,
        model_key: `no_metrics_model`,
      },
      // Model with empty metrics object
      {
        model_name: `Empty Metrics Model`,
        model_key: `empty_metrics_model`,
        metrics: {},
      },
      // Model with extreme values
      {
        model_name: `Extreme Values Model`,
        model_key: `extreme_model`,
        metrics: {
          discovery: {
            full_test_set: {
              F1: Number.MAX_VALUE,
              Accuracy: Number.MIN_VALUE,
              R2: -Infinity,
              RMSE: Infinity,
            },
          },
          phonons: {
            kappa_103: { κ_SRME: 0 },
          },
        },
      },
      // Model with all undefined values for metrics
      {
        model_name: `Undefined Metrics Model`,
        model_key: `undefined_model`,
        metrics: {
          discovery: {
            full_test_set: {
              F1: undefined,
              Accuracy: undefined,
            },
          },
          phonons: {
            kappa_103: { κ_SRME: undefined },
          },
        },
      },
    ] as unknown as ModelData[]
  }

  it(`sorts models by numeric metrics correctly with NaN handling`, () => {
    const test_models = create_test_models()

    // Test cases for different metrics and sort orders
    const test_cases = [
      {
        metric: `F1`,
        order: `desc` as const,
        expected_order: [`aaa_model`, `mmm_model`, `zzz_model`, `missing_model`],
      },
      {
        metric: `Accuracy`,
        order: `asc` as const,
        expected_order: [`zzz_model`, `aaa_model`, `mmm_model`, `missing_model`],
      },
      {
        metric: `κ_SRME`,
        order: `asc` as const,
        expected_order: [`zzz_model`, `mmm_model`, `aaa_model`, `missing_model`],
      },
    ]

    for (const { metric, order, expected_order } of test_cases) {
      const sorted_models = sort_models(test_models, metric, order)

      // Verify the order matches expected
      expected_order.forEach((model_key, idx) => {
        expect(sorted_models[idx].model_key).toBe(model_key)
      })
    }
  })

  it(`sorts models by model_name correctly`, () => {
    const test_models = create_test_models()

    // Test alphabetical sorting
    const sorted_by_name_asc = sort_models(test_models, `model_name`, `asc`)
    const sorted_by_name_desc = sort_models(test_models, `model_name`, `desc`)

    // Check ascending order
    expect(sorted_by_name_asc.map((m) => m.model_key)).toEqual([
      `aaa_model`,
      `missing_model`,
      `mmm_model`,
      `zzz_model`,
    ])

    // Check descending order
    expect(sorted_by_name_desc.map((m) => m.model_key)).toEqual([
      `zzz_model`,
      `mmm_model`,
      `missing_model`,
      `aaa_model`,
    ])
  })

  it(`handles edge cases with missing or extreme metric values`, () => {
    const edge_case_models = create_edge_case_models()
    const regular_models = create_test_models()
    const combined_models = [...regular_models, ...edge_case_models]

    // Test sorting with κ_SRME where one model has zero value
    const sorted_by_kappa = sort_models(combined_models, `κ_SRME`, `asc`)

    // Zero value should be first for asc
    expect(sorted_by_kappa[0].model_key).toBe(`extreme_model`)
  })

  it(`handles sorting with all models having missing metric`, () => {
    const all_missing_models = [
      { model_name: `Model A`, model_key: `model_a` },
      { model_name: `Model B`, model_key: `model_b` },
    ] as unknown as ModelData[]

    // Sort by a metric that none of the models have
    const sorted_models = sort_models(all_missing_models, `F1`, `desc`)

    // Order should be preserved when all models are missing the metric
    expect(sorted_models[0].model_key).toBe(`model_a`)
    expect(sorted_models[1].model_key).toBe(`model_b`)
  })

  it(`maintains original order for equivalent values`, () => {
    const models_with_same_values = [
      {
        model_name: `Model 1`,
        model_key: `model_1`,
        metrics: { discovery: { full_test_set: { F1: 0.8 } } },
      },
      {
        model_name: `Model 2`,
        model_key: `model_2`,
        metrics: { discovery: { full_test_set: { F1: 0.8 } } },
      },
      {
        model_name: `Model 3`,
        model_key: `model_3`,
        metrics: { discovery: { full_test_set: { F1: 0.8 } } },
      },
    ] as unknown as ModelData[]

    // Sort models with identical F1 values
    const sorted_models = sort_models(models_with_same_values, `F1`, `desc`)

    // Original order should be preserved
    expect(sorted_models[0].model_key).toBe(`model_1`)
    expect(sorted_models[1].model_key).toBe(`model_2`)
    expect(sorted_models[2].model_key).toBe(`model_3`)
  })
})
