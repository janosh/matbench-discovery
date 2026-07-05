import { MODELS } from '$lib'
import { get_org_logo } from '$lib/labels'
import type { ModelData } from '$lib/types'
import ModelPage from '$routes/models/[slug]/+page.svelte'
import { tick } from 'svelte'
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest'
import { mount } from '../index'

beforeEach(() =>
  vi.stubGlobal(
    `fetch`,
    vi.fn(() => Promise.resolve(new Response(`missing`, { status: 404 }))),
  ),
)
afterEach(() => vi.unstubAllGlobals())

const mirror_physics_model = MODELS.find((model) =>
  model.authors.some((author) => author.affiliation === `Mirror Physics`),
)
if (!mirror_physics_model) throw new Error(`missing Mirror Physics model`)
const test_model = mirror_physics_model

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
    expect(
      meta_info?.textContent?.includes(`Ensemble ${test_model.n_estimators} models`),
    ).toBe(test_model.n_estimators > 1)

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
      expect(author_elem.querySelector(`img.org-logo`)?.getAttribute(`src`)).toBe(
        yaml_author.affiliation ? get_org_logo(yaml_author.affiliation)?.src : undefined,
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
    expect(model_info?.textContent).toContain(test_model.targets.replaceAll(`_`, ``))
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

  it(`renders without crashing when training_set has an unknown dataset key`, () => {
    // regression: an unknown training_set key (e.g. a model-YAML typo) used to throw
    // on DATASETS[key] and crash the whole page
    const bogus_sets = [`BogusDataset`, `MPtrj`] as ModelData[`training_set`]
    const model = { ...test_model, training_set: bogus_sets }
    mount(ModelPage, { target: document.body, props: { data: { model } } })

    expect(document.querySelector(`h1`)?.textContent).toBe(test_model.model_name)
    const training_set = document.querySelector(`.training-set`)
    // known dataset still renders as link, unknown one falls back to plain text
    expect(training_set?.querySelectorAll(`a`)).toHaveLength(1)
    expect(training_set?.textContent).toContain(`BogusDataset (unknown dataset)`)
  })

  it(`lazy-mounts energy parity tab plots`, async () => {
    mount(ModelPage, { target: document.body, props: { data: { model: test_model } } })
    await tick()

    // only the default tab's plot mounts on load; toggling mounts the other for good
    expect(document.querySelectorAll(`section.energy-parity-plot`)).toHaveLength(1)
    const tab_buttons = document.querySelectorAll<HTMLButtonElement>(
      `.energy-parity-tabs button`,
    )
    tab_buttons[1].click()
    await tick()
    tab_buttons[0].click()
    await tick()
    expect(document.querySelectorAll(`section.energy-parity-plot`)).toHaveLength(2)
  })
})
