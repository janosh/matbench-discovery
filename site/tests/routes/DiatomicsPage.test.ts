import { DEFAULT_CDS_CONFIG } from '#lib/combined-scores.svelte.js'
import { ACTIVE_MODELS, has_diatomics_curves } from '#lib/models.svelte.js'
import DiatomicsPage from '#routes/benchmarks/diatomics/+page.svelte'
import { tick } from 'svelte'
import {
  DiatomicCurveCache,
  MODEL_COLORS,
} from '#routes/benchmarks/diatomics/diatomics.svelte.js'
import { describe, expect, it, vi } from 'vite-plus/test'
import {
  button_for,
  doc_query,
  expect_table_url_state,
  expect_task_page_layout,
  header_name,
  mount_with_url,
  open_filter_menu,
  patch_model,
  query_param,
  request_url,
  selected_chips,
  selected_only_toggle,
} from '../index.js'

// Real model entries with overridden names/CDS: the page resolves its model keys from
// the live MODELS, so defaults must follow these (not serialized) scores.
const fixture_scores = { a: 0.1, b: 0.8, c: 0.6, d: 0.9, e: 0.2, f: 0.3 }
const fixture_models = ACTIVE_MODELS.filter(has_diatomics_curves).slice(0, 6)
const model_keys = fixture_models.map(({ model_key }) => model_key)
const [key_a, key_b] = model_keys

const curve = {
  distances: [1],
  'homo-nuclear': {
    'H-H': { energies: [0] },
    'He-He': { energies: [0] },
    'F-F': { energies: [0] },
  },
}
const page_data = { diatomic_model_keys: model_keys, reference_names: [`PBE`, `r2SCAN`] }
const curve_url = (key: string) => `/benchmarks/diatomics/curves/${key}.json`

// serve curves/<key>.json for every key except `failing_keys` (HTTP 503)
const stub_curves = (failing_keys: string[] = []) => {
  const fetch_mock = vi.fn(async (url: RequestInfo | URL) => {
    const key = /curves\/(?<key>.+)\.json$/.exec(request_url(url))?.groups?.key
    if (!key || failing_keys.includes(key)) return new Response(`down`, { status: 503 })
    return Response.json(curve)
  })
  vi.stubGlobal(`fetch`, fetch_mock)
  return fetch_mock
}
const fetched_urls = (fetch_mock: ReturnType<typeof stub_curves>) =>
  fetch_mock.mock.calls.map(([url]) => request_url(url))

const model_select = () => doc_query(`.controls .multiselect`)

const mount_page = async (search = ``, failing_keys: string[] = []) => {
  const fetch_mock = stub_curves(failing_keys)
  for (const [idx, [suffix, combined_score]] of Object.entries(
    fixture_scores,
  ).entries()) {
    const model = fixture_models[idx]
    const { metrics } = model
    patch_model(model, {
      model_name: `Model ${suffix.toUpperCase()}`,
      metrics: { ...metrics, diatomics: { ...metrics?.diatomics, combined_score } },
    })
  }
  await mount_with_url(DiatomicsPage, `http://localhost/benchmarks/diatomics${search}`, {
    props: { data: page_data },
  })
  return fetch_mock
}

const selected_options = () =>
  doc_query(`ul[aria-label="selected options"]`, model_select())
const selected_labels = () => selected_chips(model_select())

async function select_model_option(model_name: string): Promise<void> {
  model_select().querySelector<HTMLInputElement>(`input`)?.click()
  await tick()
  const option = [
    ...document.querySelectorAll<HTMLElement>(`ul[role="listbox"] > li`),
  ].find((candidate) => candidate.textContent?.includes(model_name))
  if (!option) throw new Error(`No model option found for ${model_name}`)
  expect(option.style.color, `${model_name} option row color`).toBe(``)
  option.click()
  await tick()
}

describe(`DiatomicCurveCache`, () => {
  it(`fetches each key once, records failures and clears them on a successful retry`, async () => {
    const fetch_mock = stub_curves([`bad`])
    const cache = new DiatomicCurveCache()
    cache.request([`PBE`, `bad`, `PBE`])
    cache.request([`PBE`]) // already requested: no refetch
    await vi.waitFor(() => expect(cache.errors.get(`bad`)).toBe(`HTTP 503`))
    expect(cache.curves.get(`PBE`)).toEqual(curve)
    expect(fetched_urls(fetch_mock)).toEqual([curve_url(`PBE`), curve_url(`bad`)])

    stub_curves() // the failing key recovers
    cache.request([`bad`])
    expect(cache.errors.has(`bad`)).toBe(false) // a retry clears the stale error
    await vi.waitFor(() => expect(cache.curves.get(`bad`)).toEqual(curve))
    expect(cache.errors.has(`bad`)).toBe(false)
  })
})

describe(`Diatomics Page URL state`, () => {
  it(`defaults to the top three CDS models with curves plus DFT references`, async () => {
    const fetch_mock = await mount_page()

    expect_task_page_layout(
      [
        `Leaderboard`,
        `CDS vs twall`,
        `Test set: Diatomic curves`,
        `Diatomic Energy Curves`,
        `Methodology`,
      ],
      Object.keys(DEFAULT_CDS_CONFIG).length,
    )
    const headers = document.querySelectorAll(`section.full-bleed thead th`)
    expect(header_name([...headers].at(-1) ?? document.body)).toBe(`Reproducible`)
    expect(selected_labels()).toEqual([
      `PBE (DFT)`,
      `r2SCAN (DFT)`,
      `Model B`,
      `Model C`,
      `Model D`,
    ])
    // only the selected curves are fetched, each once
    const [key_d, key_c] = [model_keys[3], model_keys[2]]
    const default_urls = [`PBE`, `r2SCAN`, key_b, key_c, key_d].map(curve_url)
    expect(fetched_urls(fetch_mock).toSorted()).toEqual(default_urls.toSorted())
    await select_model_option(`Model A`)
    expect(selected_options()?.textContent).toContain(`Model A`)
    expect(fetched_urls(fetch_mock).toSorted()).toEqual(
      [...default_urls, curve_url(key_a)].toSorted(),
    )
    doc_query<HTMLAnchorElement>(
      `section[aria-labelledby="test-set"] a[href="#diatomic-energy-curves"]`,
    ).click()
    await tick()
    expect(selected_labels()).toEqual([`PBE (DFT)`, `r2SCAN (DFT)`])
    expect(query_param(`models`)).toBe(`PBE,r2SCAN`)
  })

  it(`restores selected curve models from the models query param and syncs edits`, async () => {
    await mount_page(`?models=unknown,${key_b},${key_b}`)

    expect(selected_labels()).toEqual([`Model B`])
    expect(query_param(`models`)).toBe(key_b)
    const selected_item = doc_query<HTMLLIElement>(`li`, selected_options())
    // palette slots follow CDS rank: Model B (0.8) is second after Model D (0.9)
    expect(selected_item.style.background).toBe(MODEL_COLORS[1])
    expect(selected_item.style.color).not.toBe(``)

    await select_model_option(`Model A`)
    expect(query_param(`models`)).toBe(`${key_a},${key_b}`)
  })

  it(`restores selected element subset from the elements query param`, async () => {
    await mount_page(`?elements=halogen`)

    expect(button_for(`All`).getAttribute(`aria-checked`)).toBe(`false`)
    expect(button_for(`Halogens`).getAttribute(`aria-checked`)).toBe(`true`)
    await vi.waitFor(() =>
      expect(document.querySelectorAll(`.diatomic-plot-title`)).toHaveLength(1),
    )
    expect(doc_query(`.diatomic-plot-title`).textContent).toContain(`F-F`)
  })

  it(`names models whose curves failed to load and still plots the rest`, async () => {
    await mount_page(`?models=PBE,${key_b}`, [key_b])
    await vi.waitFor(() =>
      expect(doc_query(`.curve-status[role="alert"]`).textContent).toContain(`Model B`),
    )
    expect(document.querySelectorAll(`.diatomic-plot-shell`).length).toBeGreaterThan(0)
  })

  it(`restores metrics-table sort and filters and syncs filter changes`, async () => {
    await mount_page(`?sort=pbe_force_mae&dir=desc&train=MPtrj&heatmap=0&selected_only=1`)

    const heatmap = await expect_table_url_state(`FMAE`, `descending`)
    const training_filter = doc_query<HTMLInputElement>(
      `[aria-label="require MPtrj"]`,
      await open_filter_menu(`Training data`),
    )
    const selected_only = selected_only_toggle()
    if (!selected_only) throw new Error(`No show-selected-only toggle`)
    expect(training_filter.checked).toBe(true)
    expect(selected_only.checked).toBe(true)

    training_filter.click()
    heatmap.click()
    selected_only.click()
    await tick()
    const params = new URL(location.href).searchParams
    for (const key of [`train`, `heatmap`, `selected_only`])
      expect(params.has(key)).toBe(false)
    expect(params.get(`sort`)).toBe(`pbe_force_mae`)
    expect(document.querySelectorAll(`tbody tr`).length).toBeGreaterThan(0)
  })

  it(`preserves metrics-table sort params when other controls update the URL`, async () => {
    await mount_page(`?models=${key_b}&sort=pbe_force_mae&dir=asc`)

    button_for(`Nonmetals`).click()
    await tick()

    const params = new URL(location.href).searchParams
    expect(params.get(`models`)).toBe(key_b)
    expect(params.get(`elements`)).toBe(`nonmetal`)
    expect(params.get(`sort`)).toBe(`pbe_force_mae`)
    // asc is non-default (default sort is CDS desc), so dir must survive in the URL
    expect(params.get(`dir`)).toBe(`asc`)
  })
})
