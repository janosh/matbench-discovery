import { goto } from '$app/navigation'
import { MODELS, score_weight_records } from '$lib/models.svelte'
import {
  calculate_cds,
  CDS_CONFIG,
  CMDS_CONFIG,
  CPS_CONFIG,
} from '$lib/combined-scores.svelte'
import { comparison } from '$lib/model-comparison.svelte'
import Layout from '$routes/+layout.svelte'
import pkg from '$site/package.json'
import { createRawSnippet, tick } from 'svelte'
import { expect, it, vi } from 'vitest'
import { doc_query, mount, mount_with_url, navigate, query_param } from '../index'

it(`restores and shares every score's weights on model, task and comparison routes`, async () => {
  const model = MODELS.find(
    ({ metrics }) =>
      metrics?.discovery?.unique_prototypes?.F1 != null &&
      metrics.md?.vdos_error != null &&
      metrics.diatomics?.run_time_sec != null,
  )
  if (!model?.metrics?.md || !model.metrics.diatomics)
    throw new Error(`No model covers the combined scores`)
  const custom_query = `cps_weights=1,0,0&cmds_weights=1,0,0,0&cds_weights=0,0,1,0`
  await mount_with_url(
    Layout,
    `http://localhost/models/${model.model_key}?${custom_query}`,
  )
  const expected_weights = {
    CPS: { F1: 1, κ_SRME: 0, RMSD: 0 },
    CMDS: { vdos_error: 1, adf_error: 0, run_time_sec: 0, pressure_error: 0 },
    CDS: { accuracy: 0, geometry: 0, speed: 1, physicality: 0 },
  }
  for (const route of [
    `/models/${model.model_key}`,
    `/benchmarks/md`,
    `/benchmarks/diatomics`,
    `/`,
  ]) {
    await navigate(`${route}?${custom_query}&compare=${model.model_key}`)
    expect(score_weight_records()).toEqual(expected_weights)
    expect(model.CPS).toBe(model.metrics.discovery?.unique_prototypes?.F1)
    expect(model.metrics.md.combined_score).toBe(
      Math.max(0, 1 - (model.metrics.md.vdos_error ?? NaN) / 100),
    )
    expect(model.metrics.diatomics.combined_score).toBe(
      calculate_cds(model.metrics.diatomics, CDS_CONFIG),
    )
    expect(query_param(`cps_weights`)).toBe(`1,0,0`)
    expect(query_param(`cmds_weights`)).toBe(`1,0,0,0`)
    expect(query_param(`cds_weights`)).toBe(`0,0,1,0`)
  }
  // Weight edits persist even when the active route has no weight-adjustment chart.
  await navigate(`/models/${model.model_key}?${custom_query}`)
  CPS_CONFIG.F1.weight = 0.12345678901234566
  CPS_CONFIG.RMSD.weight = 1 - CPS_CONFIG.F1.weight
  await tick()
  expect(query_param(`cps_weights`)).toBe(`0.12345678901234566,0,0.8765432109876543`)
  const copied_weights = score_weight_records()
  copied_weights.CMDS.vdos_error = 0
  expect(CMDS_CONFIG.vdos_error.weight).toBe(1)

  // Missing or malformed params reset their own score, rather than retaining
  // another task's in-memory weights or interpreting an unrelated score's vector.
  await navigate(`/benchmarks/diatomics?cps_weights=invalid&cmds_weights=1,0,0,0`)
  expect(query_param(`cps_weights`)).toBeNull()
  expect(query_param(`cds_weights`)).toBeNull()
  expect(query_param(`cmds_weights`)).toBe(`1,0,0,0`)
  await navigate(`/`)
  expect(location.search).toBe(``)
  comparison.open = false
  comparison.keys.clear()
})

it(`loads comparison on demand and retains its controls between openings`, async () => {
  comparison.keys.clear()
  comparison.open = false
  await mount_with_url(Layout, `http://localhost/`)
  expect(
    doc_query(`footer a[href="${pkg.repository}/blob/main/license"]`).textContent,
  ).toBe(`2022`)
  expect(document.querySelector(`dialog[aria-label="Model comparison"]`)).toBeNull()
  comparison.open_with(MODELS[0].model_key)
  const dialog = await vi.waitFor(
    () => doc_query<HTMLDialogElement>(`dialog[aria-label="Model comparison"]`),
    { timeout: 5000 },
  )
  const select = doc_query<HTMLSelectElement>(`select`, dialog)
  select.value = `added`
  select.dispatchEvent(new Event(`change`, { bubbles: true }))
  comparison.open = false
  await tick()
  comparison.open = true
  await tick()
  expect(
    doc_query<HTMLSelectElement>(`dialog[aria-label="Model comparison"] select`).value,
  ).toBe(`added`)
  comparison.open = false
  comparison.keys.clear()
})

it.each([`/api`, `/data`, `/benchmarks/diatomics`])(
  `shows the table of contents only with enough section headings on %s`,
  async (route) => {
    await mount_with_url(Layout, `http://localhost${route}`, {
      props: {
        children: createRawSnippet(() => ({
          render: () =>
            `<section><h1>Overview</h1><h2>First</h2><h2>Second</h2>${route === `/benchmarks/diatomics` ? `<h2>Third</h2>` : ``}</section>`,
        })),
      },
    })
    await vi.waitFor(() => {
      expect(document.querySelector(`.toc button`) !== null).toBe(route !== `/data`)
    })
    if (route === `/data`) return

    expect(doc_query(`.toc`).style.zIndex).toBe(`1`)
    doc_query<HTMLButtonElement>(`.toc button`).click()
    await tick()
    expect(
      [...document.querySelectorAll(`.toc nav a`)].map((link) => link.textContent),
    ).toEqual(
      route === `/api` ? [`Overview`, `First`, `Second`] : [`First`, `Second`, `Third`],
    )
    expect(doc_query(`.toc nav`).style.fontSize).toBe(`7pt`)
    expect(doc_query(`.toc-title`).style.margin).toBe(`3pt`)
  },
)

it.each([
  `/benchmarks`,
  `/benchmarks/discovery/tmi`,
  `/data/sets`,
  `/data/tmi`,
  `/models/${MODELS[0].model_key}`,
])(`navigates to %s from the command menu`, async (route) => {
  mount(Layout, { target: document.body })
  await tick()
  expect(
    [...document.querySelectorAll(`nav[data-nav] a`)].map((link) =>
      link.getAttribute(`href`),
    ),
  ).toEqual([
    `/`,
    `/benchmarks`,
    ...[`diatomics`, `discovery`, `geo-opt`, `md`, `phonons`].map(
      (task) => `/benchmarks/${task}`,
    ),
    `/models`,
    `/api`,
    `/contribute`,
    `/data/sets`,
    pkg.paper,
  ])
  expect(doc_query(`nav a[href="/benchmarks"]`).textContent).toContain(`Benchmarks`)
  expect(doc_query(`nav a[href="/data/sets"]`).textContent).toContain(`Datasets`)
  window.dispatchEvent(new KeyboardEvent(`keydown`, { key: `k`, ctrlKey: true }))
  await tick()
  const menu = document.querySelector(`dialog[aria-label="Command menu"]`)
  const option = [...(menu?.querySelectorAll(`li[role="option"]`) ?? [])].find(
    (element) => element.textContent?.trim() === route,
  )
  expect(option).toBeDefined()
  option?.dispatchEvent(new MouseEvent(`click`, { bubbles: true }))
  await tick()
  expect(goto).toHaveBeenLastCalledWith(route)
})
