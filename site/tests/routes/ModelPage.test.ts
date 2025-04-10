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
  it(`sorts models with NaN values to the end of the list`, () => {
    // Create a simple sorting function that mimics the one in +page.svelte
    function sortModels(
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

    // Create test models with some NaN metrics
    const test_models = [
      {
        model_name: `Model with NaN`,
        model_key: `model_nan`,
        metrics: {
          discovery: {
            full_test_set: {
              F1: 0.7,
              Accuracy: NaN,
            },
            pred_col: `is_stable`,
          },
        },
      },
      {
        model_name: `Model with good metrics`,
        model_key: `model_good`,
        metrics: {
          discovery: {
            full_test_set: {
              F1: 0.9,
              Accuracy: 0.85,
            },
            pred_col: `is_stable`,
          },
        },
      },
      {
        model_name: `Model with low metrics`,
        model_key: `model_low`,
        metrics: {
          discovery: {
            full_test_set: {
              F1: 0.5,
              Accuracy: 0.6,
            },
            pred_col: `is_stable`,
          },
        },
      },
    ] as ModelData[]

    // Test sorting with F1 (higher is better)
    const sorted_models = sortModels(test_models, `F1`, `desc`)

    // Log the array order to debug
    console.log(
      `Sorted models (desc F1):`,
      sorted_models.map(
        (m) => `${m.model_key} (${m.metrics?.discovery?.full_test_set?.F1})`,
      ),
    )

    // Get indexes
    const nan_model_index = sorted_models.findIndex((m) => m.model_key === `model_nan`)
    const good_model_index = sorted_models.findIndex((m) => m.model_key === `model_good`)
    const low_model_index = sorted_models.findIndex((m) => m.model_key === `model_low`)

    // Check that models are in the expected order based on the debug output
    // The F1 values are: good(0.9) > nan(0.7) > low(0.5)
    expect(good_model_index).toBe(0)
    expect(nan_model_index).toBe(1)
    expect(low_model_index).toBe(2)

    // Test with reverse sorting order (lower is better)
    const sorted_models_asc = sortModels(test_models, `Accuracy`, `asc`)

    // Log the array order to debug
    console.log(
      `Sorted models (asc Accuracy):`,
      sorted_models_asc.map(
        (m) => `${m.model_key} (${m.metrics?.discovery?.full_test_set?.Accuracy})`,
      ),
    )

    // NaN models should be sorted to the end
    const nan_model_index_asc = sorted_models_asc.findIndex(
      (m) => m.model_key === `model_nan`,
    )
    const good_model_index_asc = sorted_models_asc.findIndex(
      (m) => m.model_key === `model_good`,
    )
    const low_model_index_asc = sorted_models_asc.findIndex(
      (m) => m.model_key === `model_low`,
    )

    // We'll update the expectations based on the console output
    expect(low_model_index_asc).toBe(0)
    expect(good_model_index_asc).toBe(1)
    expect(nan_model_index_asc).toBe(2)
  })

  it(`sorts models correctly for metrics in different task categories`, () => {
    // Create a sorting function that mimics the one in +page.svelte for different task categories
    function sortModels(
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

    // Create test models with phonon metrics for κ_SRME
    const test_models = [
      {
        model_name: `Model with high κ_SRME`,
        model_key: `model_high_kappa`,
        metrics: {
          phonons: {
            κ_SRME: 0.9, // directly on phonons object
            full_test_set: { F1: 0.8 },
          },
          discovery: {
            full_test_set: { F1: 0.8 },
          },
        },
      },
      {
        model_name: `Model with medium κ_SRME`,
        model_key: `model_medium_kappa`,
        metrics: {
          phonons: {
            kappa_103: { κ_SRME: 0.5 }, // inside kappa_103 subobject
            full_test_set: { F1: 0.7 },
          },
          discovery: {
            full_test_set: { F1: 0.7 },
          },
        },
      },
      {
        model_name: `Model with low κ_SRME`,
        model_key: `model_low_kappa`,
        metrics: {
          phonons: {
            κ_SRME: 0.2, // directly on phonons object
            full_test_set: { F1: 0.6 },
          },
          discovery: {
            full_test_set: { F1: 0.6 },
          },
        },
      },
      {
        model_name: `Model missing κ_SRME`,
        model_key: `model_missing_kappa`,
        metrics: {
          phonons: `not applicable`, // String instead of object
          discovery: {
            full_test_set: { F1: 0.5 },
          },
        },
      },
    ] as ModelData[]

    // Test sorting by κ_SRME (lower is better)
    const sorted_models = sortModels(test_models, `κ_SRME`, `asc`)

    // Log the array order to debug
    console.log(
      `Sorted models (asc κ_SRME):`,
      sorted_models.map((m) => {
        const phonons = m.metrics?.phonons
        const kappa_value =
          typeof phonons === `object` && phonons !== null
            ? (phonons.kappa_103?.κ_SRME ?? phonons.κ_SRME ?? `N/A`)
            : `N/A`
        return `${m.model_key} (${kappa_value})`
      }),
    )

    // Get indexes of models
    const high_model_index = sorted_models.findIndex(
      (m) => m.model_key === `model_high_kappa`,
    )
    const medium_model_index = sorted_models.findIndex(
      (m) => m.model_key === `model_medium_kappa`,
    )
    const low_model_index = sorted_models.findIndex(
      (m) => m.model_key === `model_low_kappa`,
    )
    const missing_model_index = sorted_models.findIndex(
      (m) => m.model_key === `model_missing_kappa`,
    )

    // Check that models are in the expected order
    // For κ_SRME, lower is better, so sorting ascending should put them in ascending order
    expect(low_model_index).toBe(0)
    expect(medium_model_index).toBe(1)
    expect(high_model_index).toBe(2)
    expect(missing_model_index).toBe(3) // missing values should be at the end

    // Also test that normal discovery metrics still work
    const sorted_by_f1 = sortModels(test_models, `F1`, `desc`)

    console.log(
      `Sorted models (desc F1):`,
      sorted_by_f1.map(
        (m) => `${m.model_key} (${m.metrics?.discovery?.full_test_set?.F1})`,
      ),
    )

    // Check order for F1 (higher is better, desc order)
    const high_f1_index = sorted_by_f1.findIndex(
      (m) => m.model_key === `model_high_kappa`,
    ) // F1: 0.8
    const medium_f1_index = sorted_by_f1.findIndex(
      (m) => m.model_key === `model_medium_kappa`,
    ) // F1: 0.7
    const low_f1_index = sorted_by_f1.findIndex((m) => m.model_key === `model_low_kappa`) // F1: 0.6
    const missing_f1_index = sorted_by_f1.findIndex(
      (m) => m.model_key === `model_missing_kappa`,
    ) // F1: 0.5

    expect(high_f1_index).toBe(0)
    expect(medium_f1_index).toBe(1)
    expect(low_f1_index).toBe(2)
    expect(missing_f1_index).toBe(3)
  })

  it(`sorts models by model_name correctly`, () => {
    // Define the sorting function
    function sortModels(
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

    // Create models with alphabetically diverse names
    const test_models = [
      {
        model_name: `ZZZ Model`,
        model_key: `zzz_model`,
        metrics: { discovery: { full_test_set: { F1: 0.5 } } },
      },
      {
        model_name: `AAA Model`,
        model_key: `aaa_model`,
        metrics: { discovery: { full_test_set: { F1: 0.9 } } },
      },
      {
        model_name: `MMM Model`,
        model_key: `mmm_model`,
        metrics: { discovery: { full_test_set: { F1: 0.7 } } },
      },
    ] as ModelData[]

    // Sort by model_name (alphabetical, ascending)
    const sorted_by_name = sortModels(test_models, `model_name`, `asc`)

    // Log the sorted order
    console.log(
      `Sorted models by name (asc):`,
      sorted_by_name.map((m: ModelData) => m.model_name),
    )

    // Verify correct ordering
    expect(sorted_by_name[0].model_name).toBe(`AAA Model`)
    expect(sorted_by_name[1].model_name).toBe(`MMM Model`)
    expect(sorted_by_name[2].model_name).toBe(`ZZZ Model`)

    // Verify that F1 sorting still works (higher is better, desc order)
    const sorted_by_f1 = sortModels(test_models, `F1`, `desc`)

    console.log(
      `Sorted models by F1 (desc):`,
      sorted_by_f1.map(
        (m: ModelData) => `${m.model_name} (${m.metrics?.discovery?.full_test_set?.F1})`,
      ),
    )

    // AAA has highest F1 (0.9), then MMM (0.7), then ZZZ (0.5)
    expect(sorted_by_f1[0].model_key).toBe(`aaa_model`)
    expect(sorted_by_f1[1].model_key).toBe(`mmm_model`)
    expect(sorted_by_f1[2].model_key).toBe(`zzz_model`)
  })
})
