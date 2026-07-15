import { MODELS } from '$lib'
import { get_org_logo } from '$lib/labels'
import { RANKED_METRICS } from '$lib/rankings'
import ModelPage from '$routes/models/[slug]/+page.svelte'
import { tick } from 'svelte'
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest'
import { mount, mount_with_url } from '../index'

beforeEach(() =>
  vi.stubGlobal(
    `fetch`,
    vi.fn(() => Promise.resolve(new Response(`missing`, { status: 404 }))),
  ),
)
afterEach(() => vi.unstubAllGlobals())

const test_model = MODELS.find((model) =>
  model.authors.some((author) => author.affiliation === `Mirror Physics`),
)
if (!test_model) throw new Error(`missing Mirror Physics model`)
const test_page_data = { model: test_model, md_per_system: null }

describe(`Model Detail Page`, () => {
  it(`renders model details correctly`, () => {
    mount(ModelPage, { target: document.body, props: { data: test_page_data } })

    // Check basic model info
    expect(document.querySelector(`h1`)?.textContent).toBe(test_model.model_name)
    expect(document.body.textContent).toContain(test_model.model_version)
    expect(document.body.textContent).toContain(test_model.dates.benchmark_added)
    if (test_model.dates.paper_published)
      expect(document.body.textContent).toContain(test_model.dates.paper_published)

    // Check meta info section
    const meta_info = document.querySelector(`.meta-info`)
    expect(meta_info?.textContent).toContain(`parameters`)
    expect(
      meta_info?.textContent?.includes(`Ensemble ${test_model.n_estimators} models`),
    ).toBe((test_model.n_estimators ?? 1) > 1)
    expect(meta_info?.textContent).not.toContain(`Missing preds`)
    const discovery_detail = document.querySelector(`section.discovery-detail`)
    expect(discovery_detail?.querySelector(`h2`)?.textContent).toContain(
      `Discovery: energy and convex hull diagnostics`,
    )
    expect(
      discovery_detail?.querySelector(`.energy-parity-controls`)?.textContent,
    ).toContain(`Missing preds`)

    // Check links section
    const links = document.querySelectorAll(`.links a`)
    const expected_link_count = [
      test_model.repo,
      test_model.paper,
      test_model.docs,
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
    const expected_role =
      test_model.targets === `E` ? `Energy predictor` : `Interatomic potential`
    expect(model_info?.textContent).toContain(expected_role)
    expect(model_info?.textContent).toContain(test_model.architecture_types.join(`, `))
    expect(model_info?.textContent).toContain(test_model.targets.replaceAll(`_`, ``))
    expect(model_info?.textContent).toContain(test_model.openness)
    expect(model_info?.textContent).toContain(`Discovery Train Task`)
    expect(model_info?.textContent).toContain(`Discovery Test Task`)

    // Check training set section: one dataset link per training_sets entry
    expect(document.querySelectorAll(`.training-set a`)).toHaveLength(
      test_model.training_sets.length,
    )

    // Hyperparams section renders iff the model declares hyperparams
    const hyperparams = document.querySelector(`.hyperparams`)
    const hyperparameter_entries = Object.entries(test_model.hyperparams ?? {})
    expect(hyperparams !== null).toBe(test_model.hyperparams != null)
    for (const [key, value] of hyperparameter_entries) {
      expect(hyperparams?.textContent).toContain(key)
      expect(hyperparams?.textContent).toContain(JSON.stringify(value))
    }

    // null md_per_system page data -> no per-system MD section
    expect(document.querySelector(`section.md-per-system`)).toBeNull()
  }, 10_000)

  it(`renders leaderboard rank card with task-prefixed metric labels`, () => {
    mount(ModelPage, { target: document.body, props: { data: test_page_data } })

    const rank_links = [...document.querySelectorAll(`.rank-card a`)]
    expect(rank_links.length).toBeGreaterThan(0)
    const link_texts = rank_links.map((link) => link.textContent?.trim() ?? ``)
    // every rendered rank reads "<task-prefixed label> #<rank> /<field size>"
    for (const text of link_texts) {
      expect(text).toMatch(/#\d+ \/\d+$/)
    }
    expect(link_texts.some((text) => text.startsWith(`Discovery F1`))).toBe(true)
    // rank links lead to leaderboard pages (subset: models may lack some tasks' metrics)
    const leaderboard_hrefs = RANKED_METRICS.map((metric) => metric.rank_href)
    for (const link of rank_links) {
      expect(leaderboard_hrefs).toContain(link.getAttribute(`href`))
    }
  })

  it(`renders per-system MD breakdown from page data`, () => {
    const md_per_system = [
      {
        system: `bulkCu_1000K`,
        temperature_kelvin: 1000,
        vdos_error: 12.3,
        n_atoms: 108,
      },
      {
        system: `anthracene_293K`,
        temperature_kelvin: 293,
        vdos_error: 45.6,
        n_atoms: 72,
      },
    ]
    mount(ModelPage, {
      target: document.body,
      props: { data: { model: test_model, md_per_system } },
    })

    const section = document.querySelector(`section.md-per-system`)
    expect(section?.textContent).toContain(`per-system breakdown`)
    expect(section?.querySelectorAll(`tbody tr`)).toHaveLength(2)
    // only columns present in the rows render (no pressure/RMSE/time cols here)
    const headers = [...(section?.querySelectorAll(`th`) ?? [])].map((th) =>
      th.textContent?.replace(/\s*[↑↓]\s*$/, ``).trim(),
    )
    expect(headers).toContain(`System`)
    expect(headers).toContain(`ΔvDOS (%)`)
    expect(headers).not.toContain(`ΔP (%)`)
  })

  it(`lazy-mounts energy parity tab plots`, async () => {
    mount(ModelPage, { target: document.body, props: { data: test_page_data } })
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

  it(`restores the active energy parity tab from the energy_tab URL param`, async () => {
    await mount_with_url(
      ModelPage,
      `http://localhost/models/${test_model.model_key}?energy_tab=each`,
      { props: { data: { model: test_model } } },
    )

    const tab_buttons = document.querySelectorAll<HTMLButtonElement>(
      `.energy-parity-tabs button`,
    )
    expect(tab_buttons[0].classList.contains(`active`)).toBe(false)
    expect(tab_buttons[1].classList.contains(`active`)).toBe(true)

    // clicking back to the default tab drops the param from the URL
    tab_buttons[0].click()
    await tick()
    expect(new URL(location.href).searchParams.get(`energy_tab`)).toBeNull()
  })
})
