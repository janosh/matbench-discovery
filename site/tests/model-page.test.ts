import { MODEL_METADATA } from '$lib'
import { beforeEach, describe, expect, it } from 'vitest'
import ModelPage from '../src/routes/models/[slug]/+page.svelte'

const model_keys = MODEL_METADATA.map((m) => m.model_key)
const model_key = model_keys[0]
const test_model = MODEL_METADATA.find((m) => m.model_key === model_key)
if (!test_model) throw `missing test model`

describe(`Model Detail Page`, () => {
  beforeEach(() => {
    document.body.innerHTML = ``
  })

  it(`renders model details correctly`, async () => {
    // First verify test model exists
    new ModelPage({ target: document.body, props: { data: { model: test_model } } })

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
