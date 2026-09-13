import DATASETS from '$data/datasets.yml'
import ModelCard from '$lib/model/ModelCard.svelte'
import { ACTIVE_MODELS, MODELS } from '$lib/models.svelte'
import { ALL_METRICS } from '$lib/labels'
import { model_metric_ranks, RANKED_METRICS } from '$lib/rankings'
import { format_num } from 'matterviz/labels'
import type { ComponentProps } from 'svelte'
import { describe, expect, it } from 'vitest'
import { doc_query, mount } from '../index'

describe(`ModelCard`, () => {
  const found_model = MODELS.find((model) => model.model_key === `mace-mp-0`)
  if (!found_model) throw new Error(`Could not find mace-mp-0 model in MODELS`)
  const model = found_model

  const metrics = [ALL_METRICS.F1, ALL_METRICS.DAF, ALL_METRICS.κ_SRME]

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
      expect(header?.querySelector(`a`)?.getAttribute(`href`)).toBe(
        `/models/${model.model_key}`,
      )
      expect(header?.querySelector(`button`)).toBeNull()

      const links = document.querySelectorAll<HTMLAnchorElement>(`nav a`)
      expect(links).toHaveLength(nav_link_count(model))
      expect(links[0].href).toBe(model.repo ?? ``)
      expect(document.querySelector(`nav button`)).toBeNull()
      expect(document.body.textContent).toContain(`Added ${model.dates.benchmark_added}`)
      if (model.dates.paper_published) {
        expect(document.body.textContent).toContain(
          `Published ${model.dates.paper_published}`,
        )
      }
      expect(document.body.textContent).toContain(
        `${format_num(model.model_params, `.3~s`)} params`,
      )
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

  it.each([
    [`MPtrj`, `1,580,395 structures from 145,923 materials`],
    [`OC22`, `9,854,504 structures from 4,286 materials`],
    [`MAD-1.6`, `362,646 structures; material count not reported`],
    [`OCx24`, `unknown number of structures from 19,406 materials`],
  ])(`handles training set display for %s`, (dataset_key, expected_counts) => {
    mount_card({ model: { ...model, training_sets: [dataset_key] } })

    const dataset = DATASETS[dataset_key]
    const link = doc_query<HTMLAnchorElement>(`.metadata a[href="/data/${dataset.slug}"]`)
    expect(link.textContent).toBe(dataset_key)
    expect(link.parentElement?.textContent).toContain(`Training data:`)
    expect(link.title).toBe(`${dataset.name}: ${expected_counts}`)
  })

  describe(`Metrics Display`, () => {
    it(`displays formatted metrics with linked leaderboard ranks`, () => {
      mount_card({ metrics: RANKED_METRICS })

      const metrics_lis = document.querySelectorAll(`.metrics li`)
      expect(metrics_lis).toHaveLength(RANKED_METRICS.length)
      expect(
        [...document.querySelectorAll(`.metric-rank`)].map((link) => [
          link.textContent?.trim(),
          link.getAttribute(`href`),
        ]),
      ).toEqual(
        model_metric_ranks(model.model_key, ACTIVE_MODELS, RANKED_METRICS).map(
          ({ metric, rank, n_models }) => [`#${rank}/${n_models}`, metric.rank_href],
        ),
      )

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
      expect(document.querySelectorAll(`.metric-rank`)).toHaveLength(0)
    })
  })
})
