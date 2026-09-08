import type { ModelData } from '$lib/types'
import { MODELS } from '$lib/models.svelte'
import DiatomicsPage from '$routes/tasks/diatomics/+page.svelte'
import { tick } from 'svelte'
import { PLOT_COLORS } from 'matterviz'
import { describe, expect, it } from 'vitest'
import { doc_query, mount_with_url, sorted_header } from '../index'

const model_data: ModelData[] = (
  [
    [`a`, 0.1],
    [`b`, 0.8],
    [`c`, 0.6],
    [`d`, 0.9],
    [`e`, 1],
    [`f`, 0.95],
  ] as const
).map(([suffix, combined_score]) => ({
  ...MODELS[0],
  model_key: `model-${suffix}`,
  model_name: `Model ${suffix.toUpperCase()}`,
  metrics: { diatomics: { combined_score } },
}))

const curve = {
  distances: [1],
  'homo-nuclear': {
    'H-H': { energies: [0] },
    'He-He': { energies: [0] },
    'F-F': { energies: [0] },
  },
}
const page_data = {
  diatomic_models: model_data,
  diatomic_curves: {
    PBE: curve,
    r2SCAN: curve,
    ...Object.fromEntries(
      model_data.slice(0, 5).map(({ model_key }) => [model_key, curve]),
    ),
  },
  errors: { 'model-e': `Curve unavailable` },
  reference_names: [`PBE`, `r2SCAN`],
}

const button_for = (text: string): HTMLButtonElement => {
  const button = [...document.querySelectorAll<HTMLButtonElement>(`button`)].find(
    (candidate) => candidate.textContent?.trim() === text,
  )
  if (!button) throw new Error(`No button found for ${text}`)
  return button
}

const model_select = () => doc_query(`.controls .multiselect`)

const mount_page = async (search = ``): Promise<void> => {
  await mount_with_url(DiatomicsPage, `http://localhost/tasks/diatomics${search}`, {
    props: { data: page_data },
  })
}

const selected_options = () =>
  doc_query(`ul[aria-label="selected options"]`, model_select())
const selected_labels = () =>
  [...selected_options().querySelectorAll(`:scope > li`)].map((item) =>
    item.textContent?.trim(),
  )

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

    expect(selected_labels()).toEqual([
      `PBE (DFT)`,
      `r2SCAN (DFT)`,
      `Model B`,
      `Model C`,
      `Model D`,
    ])
    await select_model_option(`Model A`)
    expect(selected_options()?.textContent).toContain(`Model A`)
  })

  it(`restores selected curve models from the models query param`, async () => {
    await mount_page(`?models=model-b`)

    const selected = selected_options()
    expect(selected_labels()).toEqual([`Model B`])
    const selected_item = selected?.querySelector<HTMLLIElement>(`li`)
    // Model B is the second model, so it takes the second palette entry
    expect(selected_item?.style.background).toContain(PLOT_COLORS[1])
    expect(selected_item?.style.color).not.toBe(``)
  })

  it(`syncs selected curve models back to the models query param`, async () => {
    await mount_page(`?models=model-b`)

    await select_model_option(`Model A`)

    expect(new URL(location.href).searchParams.get(`models`)).toBe(`model-a,model-b`)
  })

  it(`restores selected element subset from the elements query param`, async () => {
    await mount_page(`?elements=halogen`)

    expect(button_for(`All`).getAttribute(`aria-checked`)).toBe(`false`)
    expect(button_for(`Halogens`).getAttribute(`aria-checked`)).toBe(`true`)
    expect(document.querySelectorAll(`.diatomic-plot-title`)).toHaveLength(1)
    expect(doc_query(`.diatomic-plot-title`).textContent).toContain(`F-F`)
  })

  it(`syncs selected element subset back to the elements query param`, async () => {
    await mount_page(`?models=model-b`)

    button_for(`Nonmetals`).click()
    await tick()

    const params = new URL(location.href).searchParams
    expect(params.get(`models`)).toBe(`model-b`)
    expect(params.get(`elements`)).toBe(`nonmetal`)
  })

  it(`restores metrics-table sort from sort and dir query params`, async () => {
    await mount_page(`?sort=pbe_force_mae&dir=desc`)

    expect(sorted_header()?.textContent).toContain(`PBE F MAE`)
    expect(sorted_header()?.getAttribute(`aria-sort`)).toBe(`descending`)
  })

  it(`preserves metrics-table sort params when other controls update the URL`, async () => {
    await mount_page(`?models=model-b&sort=pbe_force_mae&dir=asc`)

    button_for(`Nonmetals`).click()
    await tick()

    const params = new URL(location.href).searchParams
    expect(params.get(`models`)).toBe(`model-b`)
    expect(params.get(`elements`)).toBe(`nonmetal`)
    expect(params.get(`sort`)).toBe(`pbe_force_mae`)
    // asc is non-default (default sort is CDS desc), so dir must survive in the URL
    expect(params.get(`dir`)).toBe(`asc`)
  })
})
