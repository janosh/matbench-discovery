import { MODELS } from '$lib'
import ModelPage from '$routes/models/[slug]/+page.svelte'
import { mount } from 'svelte'
import { describe, expect, it } from 'vitest'

const model_keys = MODELS.map((m) => m.model_key)
const model_key = model_keys[0]
const test_model = MODELS.find((m) => m.model_key === model_key)
if (!test_model) throw new Error(`missing test model`)

describe(`Model Detail Page`, () => {
  it(`renders model details correctly`, () => {
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
    expect(meta_info?.textContent?.includes(`Ensemble ${test_model.n_estimators} models`)).toBe(
      test_model.n_estimators > 1,
    )

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
      expect(
        !yaml_author.affiliation ||
          author_elem.textContent?.includes(yaml_author.affiliation),
      ).toBe(true)
      expect(Boolean(author_elem.querySelector(`[href^="mailto:"]`))).toBe(
        Boolean(yaml_author.email),
      )
      expect(Boolean(author_elem.querySelector(`[href="${yaml_author.github}"]`))).toBe(
        Boolean(yaml_author.github),
      )
      expect(Boolean(author_elem.querySelector(`[href="${yaml_author.orcid}"]`))).toBe(
        Boolean(yaml_author.orcid),
      )
    }

    // Check trained by section if present
    const expected_trainers = test_model.trained_by ?? []
    const trainers = document.querySelectorAll(`.trained-by li`)
    expect(trainers).toHaveLength(expected_trainers.length)
    for (const [idx, trainer] of expected_trainers.entries()) {
      const trainer_el = trainers[idx]
      expect(trainer_el.textContent).toContain(trainer.name)
      expect(
        !trainer.affiliation || trainer_el.textContent?.includes(trainer.affiliation),
      ).toBe(true)
    }

    // Check model info section
    const model_info = document.querySelector(`.model-info`)
    expect(model_info?.textContent).toContain(test_model.model_type)
    expect(model_info?.textContent).toContain(test_model.targets)
    expect(model_info?.textContent).toContain(test_model.openness)

    // Check training set section if present
    const training_set = document.querySelector(`.training-set`)
    const train_set_links = training_set?.querySelectorAll(`a`)
    const expected_train_set_links = test_model.training_set
      ? Array.isArray(test_model.training_set)
        ? test_model.training_set.length
        : 1
      : 0
    expect(training_set !== null || !test_model.training_set).toBe(true)
    expect(train_set_links?.length ?? 0).toBe(expected_train_set_links)

    // Check hyperparameters section if present
    const hyperparams = document.querySelector(`.hyperparams`)
    const hyperparam_entries = Object.entries(test_model.hyperparams ?? {})
    expect(hyperparams !== null || hyperparam_entries.length === 0).toBe(true)
    for (const [key, value] of hyperparam_entries) {
      expect(hyperparams?.textContent).toContain(key)
      expect(hyperparams?.textContent).toContain(JSON.stringify(value))
    }
  }, 10_000)
})
