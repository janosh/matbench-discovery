import { ACTIVE_MODELS, MODELS, update_models_cps } from '$lib/models.svelte'
import {
  CDS_CONFIG,
  CMDS_CONFIG,
  CPS_CONFIG,
  update_models_cds,
  update_models_cmds,
} from '$lib/combined-scores.svelte'
import { parse_dependency_spec } from '$lib/environment'
import { get_org_logo } from '$lib/labels'
import { model_metric_ranks, RANKED_METRICS } from '$lib/rankings'
import {
  discovery_task_tooltips,
  openness_tooltips,
  targets_tooltips,
} from '$lib/metrics'
import type { ModelData } from '$lib/types'
import * as predictions from '$lib/server/predictions'
import DATASETS from '$data/datasets.yml'
import ModelPage from '$routes/models/[slug]/+page.svelte'
import { load } from '$routes/models/[slug]/+page.server'
import { format_num } from 'matterviz/labels'
import { type ComponentProps, tick } from 'svelte'
import { describe, expect, it, onTestFinished, vi } from 'vitest'
import { doc_query, mount, mount_with_url } from '../index'

const test_model = MODELS.find((model) =>
  model.authors.some((author) => author.affiliation === `Mirror Physics`),
)
if (!test_model) throw new Error(`missing Mirror Physics model`)
type PageData = ComponentProps<typeof ModelPage>[`data`]
const fixture_data = (
  model: ModelData,
  md_per_system: PageData[`md_per_system`] = null,
): PageData => {
  const canonical = MODELS.find((entry) => entry.model_key === model.model_key)
  if (!canonical) throw new Error(`Missing fixture model ${model.model_key}`)
  const original = { ...canonical }
  Object.assign(canonical, model)
  onTestFinished(() => {
    for (const key of Object.keys(canonical))
      if (!(key in original)) Reflect.deleteProperty(canonical, key)
    Object.assign(canonical, original)
  })
  return { model_key: model.model_key, md_per_system }
}
const mount_page = (model = test_model) =>
  mount(ModelPage, { target: document.body, props: { data: fixture_data(model) } })
const locator_model = MODELS.find((model) => model.model_key === `tace-oam-l`)
if (!locator_model) throw new Error(`missing TACE-OAM-L model`)
const locator_dependency = locator_model.environment.dependencies[0]
if (!locator_dependency) throw new Error(`missing TACE-OAM-L dependency`)
const locator_detail = parse_dependency_spec(locator_dependency).detail

describe(`Model Detail Page`, () => {
  it(`keeps serialized page data tied to live CPS, CMDS and CDS values and ranks`, async () => {
    const score_model = MODELS.find(
      (entry) =>
        Number.isFinite(entry.CPS) &&
        Number.isFinite(entry.metrics?.md?.combined_score) &&
        Number.isFinite(entry.metrics?.diatomics?.combined_score),
    )
    if (!score_model) throw new Error(`Missing model with all combined scores`)
    vi.spyOn(predictions, `read_md_per_system`).mockResolvedValue(null)
    // Exercise the serialized server/client boundary instead of sharing object references.
    const serialized = JSON.stringify(
      await load({ params: { slug: score_model.model_key } } as Parameters<
        typeof load
      >[0]),
    )
    const data: PageData = JSON.parse(serialized)
    expect(data).toEqual({ model_key: score_model.model_key, md_per_system: null })
    mount(ModelPage, { target: document.body, props: { data } })
    await tick()
    const initial_values = [
      score_model.CPS,
      score_model.metrics?.md?.combined_score,
      score_model.metrics?.diatomics?.combined_score,
    ]
    for (const [config, focus, update] of [
      [CPS_CONFIG, `F1`, () => update_models_cps(MODELS, CPS_CONFIG)],
      [CMDS_CONFIG, `vdos_error`, () => update_models_cmds(MODELS, CMDS_CONFIG)],
      [CDS_CONFIG, `speed`, () => update_models_cds(MODELS, CDS_CONFIG)],
    ] as const) {
      const entries = Object.entries(config)
      const weights = entries.map(([, entry]) => entry.weight)
      onTestFinished(() => {
        entries.forEach(([, entry], idx) => (entry.weight = weights[idx]))
        update()
      })
      for (const [key, entry] of entries) entry.weight = Number(key === focus)
      update()
    }
    await tick()
    const ranks = model_metric_ranks(score_model.model_key, ACTIVE_MODELS, RANKED_METRICS)
    const cps = ranks.find(({ metric }) => metric.key === `CPS`)
    if (!cps) throw new Error(`Missing weighted CPS rank`)
    expect(cps.value).not.toBe(initial_values[0])
    expect(doc_query(`.overview a[href="/"]`).textContent).toBe(
      `CPS ${format_num(cps.value, `.3`)}`,
    )
    for (const [idx, task] of [`md`, `diatomics`].entries()) {
      const rank = ranks.find(({ metric }) => metric.rank_href === `/benchmarks/${task}`)
      if (!rank) throw new Error(`Missing weighted ${task} rank`)
      expect(rank.value).not.toBe(initial_values[idx + 1])
      const card = [...document.querySelectorAll(`.rank-card a`)].find((link) =>
        link.textContent?.includes(task === `md` ? `MD CMDS` : `Diatomics CDS`),
      )
      expect(card?.querySelector(`strong`)?.textContent).toBe(
        format_num(rank.value, rank.metric.format ?? `.3`),
      )
      expect(card?.querySelector(`small`)?.textContent).toBe(
        `#${rank.rank}/${rank.n_models}`,
      )
    }
  })

  it(`renders model details correctly`, async () => {
    const pypi = `https://pypi.org/project/test-model`
    const hyperparams = {
      training: { learning_rate: 0.001 },
      upstream_config: { layers: [32, 64], enabled: true, label: `<model> & "values"` },
    }
    const notes = {
      html: {
        description: `<h3>Model notes</h3><p>Model description</p><ul><li>Run notes</li></ul>`,
      },
    }
    mount_page({ ...test_model, hyperparams, pypi, notes })
    expect(document.querySelector(`.discovery-detail`)).toBeNull()
    await tick()
    doc_query<HTMLAnchorElement>(
      `.rank-card a[href="?diagnostic=discovery#diagnostics"]`,
    ).click()
    await tick()

    expect(document.querySelector(`h1`)?.textContent).toBe(test_model.model_name)
    expect(document.body.textContent).toContain(test_model.model_version)
    expect(document.body.textContent).toContain(test_model.dates.benchmark_added)
    expect(document.querySelector(`.notes p`)?.textContent).toBe(`Model description`)
    expect(document.querySelector(`.notes li`)?.textContent).toBe(`Run notes`)
    expect(document.querySelector(`.notes p p, .notes p ul`)).toBeNull()
    const note_heading = doc_query(`.notes h3`)
    expect(note_heading.textContent).toBe(`Model notes`)
    expect(getComputedStyle(note_heading).textAlign).not.toBe(`center`)
    expect(getComputedStyle(doc_query(`.discovery-detail h3`)).textAlign).toBe(`center`)
    if (test_model.dates.paper_published)
      expect(document.body.textContent).toContain(test_model.dates.paper_published)

    const meta_info = document.querySelector(`.meta-info`)
    expect(meta_info?.textContent).toContain(`parameters`)
    expect(meta_info?.querySelector(`code`)?.textContent).toContain(
      `uv pip install test-model`,
    )
    expect(
      meta_info?.textContent?.includes(`Ensemble ${test_model.n_estimators} models`),
    ).toBe((test_model.n_estimators ?? 1) > 1)
    expect(meta_info?.textContent).not.toContain(`Missing preds`)
    const discovery_detail = document.querySelector(`section.discovery-detail`)
    expect(
      discovery_detail?.querySelector(`h2`)?.textContent?.replaceAll(/\s+/g, ` `),
    ).toContain(`Discovery: energy and convex hull diagnostics`)
    expect(doc_query(`.energy-parity-controls .missing-preds`).textContent).toContain(
      `Missing preds: ${test_model.metrics?.discovery?.full_test_set?.missing_preds}`,
    )

    const links = document.querySelectorAll(`.links a`)
    const expected_link_count = [
      test_model.repo,
      test_model.paper,
      test_model.docs,
      test_model.doi,
      pypi,
    ].filter(Boolean).length
    expect(links.length).toBeGreaterThanOrEqual(expected_link_count)
    const rank_layout = getComputedStyle(doc_query(`.rank-card`))
    expect(rank_layout.rowGap).toBe(rank_layout.columnGap)
    for (const control of document.querySelectorAll(
      `.links :is(a, summary, button), .rank-card a`,
    )) {
      const style = getComputedStyle(control)
      expect([``, `auto`, `0px`]).toContain(style.minHeight)
      expect(style.paddingTop).toBe(style.paddingBottom)
    }

    const authors = document.querySelectorAll(`.authors li`)
    expect(authors).toHaveLength(test_model.authors.length)
    for (const [idx, yaml_author] of test_model.authors.entries()) {
      const author_elem = authors[idx]
      expect(author_elem.textContent).toContain(yaml_author.name)
      expect(
        !yaml_author.affiliation ||
          author_elem.textContent?.includes(yaml_author.affiliation),
      ).toBe(true)
      const email_link = author_elem.querySelector<HTMLElement>(`[href^="mailto:"]`)
      expect(Boolean(email_link)).toBe(Boolean(yaml_author.email))
      if (email_link) {
        expect(email_link.style.fontSize).toBe(`1.1em`)
        expect(email_link.style.paddingRight).toBe(`0.2em`)
      }
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

    const model_info = doc_query(`.model-info`)
    const expected_role =
      test_model.targets === `E` ? `Energy predictor` : `Interatomic potential`
    expect(doc_query(`.overview`).textContent).toContain(expected_role)
    expect(model_info.querySelector(`ul`)).toBeNull()
    expect(
      [...model_info.querySelectorAll(`dl`)].map((list) =>
        list.getAttribute(`aria-label`),
      ),
    ).toEqual([`Specifications`, `Discovery protocol`])
    const info_rows = [...model_info.querySelectorAll(`dl > div`)].map((row) => [
      row.querySelector(`dt`)?.textContent,
      row.querySelector(`dd`)?.textContent?.trim().replaceAll(/\s+/g, ` `),
    ])
    expect(info_rows).toEqual([
      [`Architecture`, test_model.architecture_types.join(`, `)],
      [
        `Outputs`,
        `${targets_tooltips[test_model.targets]} ${test_model.targets.replaceAll(`_`, ``)}`,
      ],
      [`Openness`, `${openness_tooltips[test_model.openness]} ${test_model.openness}`],
      [
        `Training`,
        `${discovery_task_tooltips[test_model.train_task]} ${test_model.train_task}`,
      ],
      [
        `Evaluation`,
        `${discovery_task_tooltips[test_model.test_task]} ${test_model.test_task}`,
      ],
    ])
    const detail_position = doc_query(`.overview`).compareDocumentPosition(
      doc_query(`.discovery-detail`),
    )
    expect(Boolean(detail_position & Node.DOCUMENT_POSITION_FOLLOWING)).toBe(true)
    expect(
      [...document.querySelectorAll(`.licenses dt`)].map((item) => item.textContent),
    ).toEqual([`Code license`, `Checkpoint license`])
    expect(
      [...document.querySelectorAll(`.licenses dd`)].map((item) =>
        item.textContent?.trim(),
      ),
    ).toEqual([test_model.license.code, test_model.license.checkpoint])

    // one dataset link per training_sets entry
    expect(document.querySelectorAll(`.training-set a`)).toHaveLength(
      test_model.training_sets.length,
    )
    for (const dataset_key of test_model.training_sets)
      expect(doc_query(`.training-set`).textContent).toContain(
        DATASETS[dataset_key].license,
      )

    const hyperparams_tree = doc_query(`.hyperparams .json-tree`)
    expect(
      hyperparams_tree.querySelector(`.json-value.number`)?.textContent?.trim(),
    ).toBe(`0.001`)
    expect(
      hyperparams_tree.querySelector(`.json-value.boolean`)?.textContent?.trim(),
    ).toBe(`true`)
    expect(hyperparams_tree.querySelector(`.json-value.string`)?.textContent).toContain(
      hyperparams.upstream_config.label,
    )
    expect(hyperparams_tree.querySelector(`model`)).toBeNull()
    expect(
      hyperparams_tree.querySelector(
        `.json-tree-header, .collapse-level-btn, .size-hint`,
      ),
    ).toBeNull()
    const training_node = doc_query(`[data-path="training"]`, hyperparams_tree)
    training_node.click()
    await tick()
    expect(hyperparams_tree.querySelector(`.path-breadcrumb`)).toBeNull()
    doc_query<HTMLButtonElement>(`.collapse-toggle`, training_node).click()
    await tick()
    expect(training_node.getAttribute(`aria-expanded`)).toBe(`false`)

    // null md_per_system page data -> no per-system MD section
    expect(document.querySelector(`section.md-per-system`)).toBeNull()
  }, 10_000)

  it.each([
    `https://pypi.org/project/test-model/`,
    `https://pypi.org/project/test-model/1.2.3/?source=docs#files`,
  ])(`builds the install command from %s`, (pypi) => {
    mount_page({ ...test_model, pypi })
    expect(document.querySelector(`.meta-info code`)?.textContent?.trim()).toBe(
      `uv pip install test-model`,
    )
  })

  it.each([true, false])(
    `shows shared-runner instructions only with declared phonon settings: %s`,
    (configured) => {
      mount_page({
        ...test_model,
        hyperparams: configured
          ? { evaluation: { kappa: { protocol: `phonondb-v1` } } }
          : undefined,
      })
      const command = document.querySelector(`.runner-command code`)
      expect(Boolean(command)).toBe(configured)
      if (command)
        expect(command.textContent).toBe(
          `uv run models/run_kappa.py --model ${test_model.model_key} --print-cmd --dry-run`,
        )
      expect(doc_query(`.run-model a[href="#dependencies"]`)).not.toBeNull()
    },
  )

  it(`preserves both ends of overflowing dependency details`, () => {
    mount_page({ ...locator_model, hyperparams: undefined })

    expect(document.querySelector(`.hyperparams`)).toBeNull()
    const link = document.querySelector<HTMLAnchorElement>(`.deps .dependency-detail`)
    if (!link) throw new Error(`missing dependency detail link`)
    const [leading, trailing] = link.querySelectorAll(`span`)
    if (!leading || !trailing) throw new Error(`missing dependency detail spans`)
    // happy-dom doesn't reflect aria-* attributes to element properties (ariaLabel)
    expect(link.getAttribute(`aria-label`)).toBe(locator_detail)
    expect(link.title).toBe(locator_detail)
    expect(leading.textContent).toBe(locator_detail.slice(0, -10))
    expect(trailing.textContent).toBe(locator_detail.slice(-10))
    expect(getComputedStyle(leading).textOverflow).toBe(`ellipsis`)
  })

  it(`navigates available task diagnostics and keeps unavailable tasks linked to their leaderboard`, async () => {
    await mount_with_url(ModelPage, `http://localhost/models/${test_model.model_key}`, {
      props: {
        data: fixture_data({
          ...test_model,
          metrics: { ...test_model.metrics, md: undefined, diatomics: undefined },
        }),
      },
    })

    const rank_links = [...document.querySelectorAll(`.rank-card a`)]
    expect(rank_links).toHaveLength(RANKED_METRICS.length - 1)
    const link_texts = rank_links.map((link) => link.textContent?.trim() ?? ``)
    expect(link_texts.some((text) => text.startsWith(`Discovery F1`))).toBe(true)
    expect(document.querySelector(`.energy-parity-plot`)).toBeNull()
    expect(rank_links.some((link) => link.textContent?.includes(`No results`))).toBe(true)
    rank_links[0]?.dispatchEvent(new MouseEvent(`mouseenter`))
    await vi.waitFor(() =>
      expect(document.querySelector(`.popover`)?.textContent?.trim()).toMatch(
        /^Ranked \d+ of \d+ active models with this metric/,
      ),
    )
    expect(document.querySelector(`.popover br`)).not.toBeNull()
    for (const [idx, link] of rank_links.entries()) {
      if (link.textContent?.includes(`No results`))
        expect(link.getAttribute(`href`)).toMatch(/^\/benchmarks\//)
      else
        expect(link.getAttribute(`href`)).toBe(
          `?diagnostic=${[`discovery`, `geo-opt`, `phonons`][idx]}#diagnostics`,
        )
    }
    rank_links[0].dispatchEvent(new MouseEvent(`click`, { bubbles: true, ctrlKey: true }))
    await tick()
    expect(document.querySelector(`.discovery-detail`)).toBeNull()
    rank_links[0].dispatchEvent(
      new MouseEvent(`click`, { bubbles: true, cancelable: true }),
    )
    await tick()
    expect(new URL(location.href).searchParams.get(`diagnostic`)).toBe(`discovery`)
    expect(document.querySelector(`.discovery-detail`)).not.toBeNull()
    expect(rank_links[0].getAttribute(`aria-current`)).toBe(`true`)
    rank_links[1].dispatchEvent(new MouseEvent(`click`, { bubbles: true }))
    await tick()
    expect(new URL(location.href).searchParams.get(`diagnostic`)).toBe(`geo-opt`)
    expect(document.querySelector(`.discovery-detail`)).toBeNull()
    expect(doc_query(`.diagnostic-context a`).getAttribute(`href`)).toBe(
      `/benchmarks/geo-opt?models=${test_model.model_key}`,
    )
    rank_links[2].dispatchEvent(new MouseEvent(`click`, { bubbles: true }))
    await tick()
    expect(doc_query(`.diagnostic-context a`).getAttribute(`href`)).toBe(
      `/benchmarks/phonons?model=${test_model.model_key}`,
    )
  })

  it(`renders per-system MD breakdown from page data`, async () => {
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
    await mount_with_url(
      ModelPage,
      `http://localhost/models/${test_model.model_key}?diagnostic=md`,
      {
        props: {
          data: fixture_data(
            {
              ...test_model,
              metrics: {
                ...test_model.metrics,
                md: {
                  pred_file: { name: `md.csv.gz`, url: `https://example.com/md.csv.gz` },
                },
              },
            },
            md_per_system,
          ),
        },
      },
    )

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

  it(`lazy-mounts energy parity plots and observes their size`, async () => {
    const observe = vi.spyOn(ResizeObserver.prototype, `observe`)
    mount_page()
    await tick()
    doc_query<HTMLAnchorElement>(
      `.rank-card a[href="?diagnostic=discovery#diagnostics"]`,
    ).click()
    await tick()

    // only the default tab's plot mounts on load; toggling mounts the other for good
    expect(document.querySelectorAll(`section.energy-parity-plot`)).toHaveLength(1)
    const tab_buttons = document.querySelectorAll<HTMLButtonElement>(
      `.energy-parity-tabs button`,
    )
    // Only loading tabs need spinner space; inactive/finished tabs fit their labels.
    expect(tab_buttons[0].querySelector(`.circle-spinner`)).not.toBeNull()
    expect(tab_buttons[1].querySelector(`.circle-spinner`)).toBeNull()
    tab_buttons[1].click()
    await tick()
    expect(tab_buttons[0].querySelector(`.circle-spinner`)).toBeNull()
    expect(tab_buttons[1].querySelector(`.circle-spinner`)).not.toBeNull()
    tab_buttons[0].click()
    await tick()
    expect(document.querySelectorAll(`section.energy-parity-plot`)).toHaveLength(2)
    expect(observe).toHaveBeenCalledWith(
      document.querySelector(`section.energy-parity-plot`),
    )
    await vi.waitFor(() => {
      expect(document.querySelector(`.energy-parity-tabs .circle-spinner`)).toBeNull()
    })
  })

  it(`restores the active energy parity tab from the energy_tab URL param`, async () => {
    await mount_with_url(
      ModelPage,
      `http://localhost/models/${test_model.model_key}?diagnostic=discovery&energy_tab=each`,
      { props: { data: { model_key: test_model.model_key, md_per_system: null } } },
    )

    expect(
      document.querySelector(`[aria-label="Energy parity diagnostics"]`),
    ).not.toBeNull()
    const tab_buttons = document.querySelectorAll<HTMLButtonElement>(
      `.energy-parity-tabs button`,
    )
    expect(tab_buttons[0].getAttribute(`aria-checked`)).toBe(`false`)
    expect(tab_buttons[1].getAttribute(`aria-checked`)).toBe(`true`)

    // clicking back to the default tab drops the param from the URL
    tab_buttons[0].click()
    await tick()
    expect(new URL(location.href).searchParams.get(`energy_tab`)).toBeNull()
  })
})
