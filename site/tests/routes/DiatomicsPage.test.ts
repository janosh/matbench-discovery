import { DEFAULT_CDS_CONFIG } from '#lib/combined-scores.svelte.js'
import { ACTIVE_MODELS, has_diatomics_curves } from '#lib/models.svelte.js'
import DiatomicsPage from '#routes/benchmarks/diatomics/+page.svelte'
import { tick } from 'svelte'
import { MODEL_COLORS } from '#routes/benchmarks/diatomics/diatomics.svelte.js'
import { describe, expect, it } from 'vite-plus/test'
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
const page_data = {
  diatomic_model_keys: model_keys,
  diatomic_curves: {
    PBE: curve,
    r2SCAN: curve,
    ...Object.fromEntries(model_keys.map((key) => [key, curve])),
  },
  reference_names: [`PBE`, `r2SCAN`],
}

const model_select = () => doc_query(`.controls .multiselect`)

const mount_page = async (search = ``): Promise<void> => {
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

describe(`Diatomics Page URL state`, () => {
  it(`defaults to the top three CDS models with curves plus DFT references`, async () => {
    await mount_page()

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
    await select_model_option(`Model A`)
    expect(selected_options()?.textContent).toContain(`Model A`)
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
    expect(document.querySelectorAll(`.diatomic-plot-title`)).toHaveLength(1)
    expect(doc_query(`.diatomic-plot-title`).textContent).toContain(`F-F`)
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
