import { DATASETS, ModelCard, MODELS } from '$lib'
import { ALL_METRICS } from '$lib/labels'
import { format_num } from 'matterviz'
import type { ComponentProps } from 'svelte'
import { describe, expect, it } from 'vitest'
import { mount } from '../index'

describe(`ModelCard`, () => {
  // Get a real model from MODELS
  const found_model = MODELS.find((model) => model.model_key === `mace-mp-0`)
  if (!found_model) throw new Error(`Could not find mace-mp-0 model in MODELS`)
  const model = found_model

  const metrics = ([`F1`, `DAF`, `κ_SRME`] as const).map((metric_key) => ({
    ...ALL_METRICS[metric_key],
    key: metric_key,
  }))

  // Mount ModelCard with the shared model/metrics, overriding props per test
  const mount_card = (overrides: Partial<ComponentProps<typeof ModelCard>> = {}) =>
    mount(ModelCard, {
      target: document.body,
      props: { model, metrics, sort_by: `F1`, ...overrides },
    })

  // nav renders one link per http(s) URL among repo/paper/docs/checkpoint_url,
  // plus the always-present Files link built from pkg.repository
  const nav_link_count = (mdl: typeof model) =>
    [mdl.repo, mdl.paper, mdl.docs, mdl.checkpoint_url].filter((href) =>
      href?.startsWith(`http`),
    ).length + 1

  describe(`Basic Rendering`, () => {
    it(`renders model header and basic info`, () => {
      mount_card()

      const header = document.querySelector(`h2`)
      expect(header?.textContent).toContain(`MACE`)

      const links = document.querySelectorAll<HTMLAnchorElement>(`nav a`)
      expect(links).toHaveLength(nav_link_count(model))
      expect(links[0].href).toBe(model.repo ?? ``)
      expect(document.body.textContent).toContain(`Added ${model.dates.benchmark_added}`)
      expect(
        !model.dates.paper_published ||
          document.body.textContent?.includes(`Published ${model.dates.paper_published}`),
      ).toBe(true)
      expect(document.body.textContent).toContain(
        `${format_num(model.model_params, `.3~s`)} params`,
      )

      // show_details defaults to false: no detail sections rendered
      expect(
        document.querySelectorAll(`section:not(.metrics):not(.metadata) h3`),
      ).toHaveLength(0)
    })

    it(`handles missing optional fields gracefully`, () => {
      const stripped = {
        ...model,
        dates: { ...model.dates, paper_published: null },
        paper: null,
        docs: null,
      }
      mount_card({ model: stripped })

      expect(document.querySelectorAll(`nav a`)).toHaveLength(nav_link_count(stripped))
      expect(document.body.textContent).not.toContain(`Published`)
    })
  })

  it(`handles training set display`, () => {
    mount_card()

    // Look for span containing "Training set" text
    const training_set = [...document.querySelectorAll(`section.metadata span`)].find(
      (span) => span.textContent?.includes(`Training data`),
    )
    expect(training_set?.textContent).toContain(`Training data:`)

    // Test actual training set data
    const training_set_links = training_set?.querySelectorAll(`a`)
    const dataset_key = model.training_sets[0]
    const dataset = DATASETS[dataset_key]

    // Check that we're linking to our internal data page
    expect(training_set_links?.[0]?.href).toContain(`/data/${dataset.slug}`)

    // Check that structure count is shown in tooltip
    const formatted_structures = format_num(dataset.n_structures)
    expect(training_set_links?.[0]?.title).toContain(`${formatted_structures} structures`)
  })

  describe(`Metrics Display`, () => {
    it(`displays metrics with correct formatting`, () => {
      mount_card()

      const metrics_lis = document.querySelectorAll(`.metrics li`)
      expect(metrics_lis).toHaveLength(metrics.length)

      const f1_metric = [...metrics_lis].find((item) => item.textContent?.includes(`F1`))
      const f1_value = model.metrics?.discovery?.unique_prototypes?.F1
      expect(f1_metric?.querySelector(`strong`)?.textContent?.trim()).toBe(
        f1_value?.toString(),
      )
      expect(f1_metric?.classList.contains(`active`)).toBe(true)

      const kappa_metric = [...metrics_lis].find((item) =>
        item.textContent?.includes(`κ`),
      )
      const kappa_value = Number(model.metrics?.phonons?.kappa_103?.κ_SRME) || 0
      const displayed_kappa = kappa_metric?.querySelector(`strong`)?.textContent?.trim()
      // must not render blank (Number(``) would coerce to 0 and mask that case)
      expect(displayed_kappa).toMatch(/\d/)
      // The displayed value may be rounded differently
      expect(Number(displayed_kappa)).toBeCloseTo(kappa_value, 2)
    })

    it(`handles missing metrics`, () => {
      const model_without_metrics = { ...model, metrics: undefined }

      mount_card({ model: model_without_metrics })

      const metrics_li_strong = document.querySelectorAll(`.metrics li strong`)[0]
      expect(metrics_li_strong.textContent?.trim()).toBe(`n/a`)
    })
  })

  describe(`Expandable Details`, () => {
    it(`displays authors and package versions correctly`, () => {
      mount_card({ show_details: true })

      // Check author info within the list item
      const author_li = document.querySelector(`section:first-child ul li`)
      expect(author_li?.textContent?.trim()).toContain(model.authors[0].name)

      // Check package versions
      const packages = [...document.querySelectorAll(`section:nth-child(2) li`)]
      const expected = model.environment.dependencies
      expect(packages).toHaveLength(expected.length)
      const first_dep = expected[0] ?? ``
      expect(packages[0]?.textContent?.trim() ?? ``).toContain(
        first_dep.split(/[=<>!~@]/)[0]?.trim() ?? ``,
      )
    })
  })
})
