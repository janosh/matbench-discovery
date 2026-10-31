import DATASETS from '$data/datasets.yml'
import ModelCard from '#lib/model/ModelCard.svelte'
import { ACTIVE_MODELS, MODELS } from '#lib/models.svelte.js'
import { ALL_METRICS, MD_METRICS } from '#lib/labels.js'
import { get_nested_number, label_data_path } from '#lib/metrics.js'
import { model_metric_ranks, RANKED_METRICS } from '#lib/rankings.js'
import { format_num } from 'matterviz/labels'
import { strip_html } from 'svelte-widgets/utils'
import type { ComponentProps } from 'svelte'
import { describe, expect, it, vi } from 'vite-plus/test'
import { doc_query, mount } from '../index.js'

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
      expect(document.body.textContent).toContain(
        `Published ${model.dates.paper_published}`,
      )
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

      // values use each metric's own format (and unit); rank links are named by the
      // metric's plain-text label, not its internal key
      for (const [idx, metric] of RANKED_METRICS.entries()) {
        const value = get_nested_number(model, label_data_path(metric))
        expect(metrics_lis[idx].querySelector(`strong`)?.textContent?.trim()).toBe(
          `${format_num(value ?? NaN, metric.format)} ${metric.unit ?? ``}`.trim(),
        )
      }
      expect(
        [...document.querySelectorAll(`.metric-rank`)].map((link) =>
          link.getAttribute(`aria-label`),
        ),
      ).toEqual(
        model_metric_ranks(model.model_key, ACTIVE_MODELS, RANKED_METRICS).map(
          ({ metric, rank, n_models }) =>
            `${strip_html(metric.label)}: rank ${rank} of ${n_models}`,
        ),
      )
      const f1_metric = [...metrics_lis].find((item) => item.textContent?.includes(`F1`))
      expect(f1_metric?.classList.contains(`active`)).toBe(true)
    })

    it.each([undefined, NaN, Infinity])(
      `explains unavailable metrics (%s)`,
      async (value) => {
        const kappa_metrics = model.metrics?.phonons?.kappa_103
        if (!kappa_metrics)
          throw new Error(`Missing phonon metrics for ${model.model_key}`)
        mount_card({
          model: {
            ...model,
            license: { ...model.license, checkpoint: `unreleased` },
            metrics:
              value === undefined
                ? {}
                : { phonons: { kappa_103: { ...kappa_metrics, κ_SRME: value } } },
          },
          metrics: [ALL_METRICS.κ_SRME],
        })

        const metrics_li_strong = document.querySelectorAll(`.metrics li strong`)[0]
        expect(metrics_li_strong.textContent?.trim()).toBe(`n/a`)
        expect(document.querySelectorAll(`.metric-rank`)).toHaveLength(0)
        doc_query(`.metrics [role="button"]`).focus()
        await vi.waitFor(() => {
          const text = doc_query(`.popover`).textContent
          expect(text).toContain(`Model weights are not publicly available.`)
          expect(text).not.toContain(`Contributions welcome`)
        })
      },
    )

    it(`names the Speed subscore GPU in the CMDS popover`, async () => {
      const cost = { source: `runner`, hardware: `H200`, run_time_sec: 2e4 } as const
      const md = { cost, combined_score: 0.6 }
      mount_card({
        model: { ...model, metrics: { md } },
        metrics: [MD_METRICS.md_combined_score],
      })
      expect(doc_query(`.metrics li strong`).textContent?.trim()).toBe(`0.600`)
      doc_query(`.metrics [role="button"]`).focus()
      await vi.waitFor(() =>
        expect(doc_query(`.popover`).textContent).toContain(
          `Speed subscore timed on H200 (benchmark runner)`,
        ),
      )
    })
  })
})
