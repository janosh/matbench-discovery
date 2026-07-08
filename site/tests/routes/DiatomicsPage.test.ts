import type { ModelData } from '$lib/types'
import DiatomicsPage from '$routes/tasks/diatomics/+page.svelte'
import { tick } from 'svelte'
import { describe, expect, it } from 'vitest'
import { mount_with_url, sorted_header } from '../index'

const model_data = [
  { model_key: `model-a`, model_name: `Model A` },
  { model_key: `model-b`, model_name: `Model B` },
] as ModelData[]

const curve = { distances: [1], 'homo-nuclear': { 'H-H': { energies: [0] } } }
const page_data = {
  diatomic_models: model_data,
  diatomic_curves: { PBE: curve, r2SCAN: curve, 'Model A': curve, 'Model B': curve },
  errors: {},
  reference_names: [`PBE`, `r2SCAN`],
}

const button_for = (text: string): HTMLButtonElement => {
  const button = [...document.querySelectorAll<HTMLButtonElement>(`button`)].find(
    (candidate) => candidate.textContent?.trim() === text,
  )
  if (!button) throw new Error(`No button found for ${text}`)
  return button
}

const model_select = (): HTMLElement => {
  // scope to the curve-controls div: the Pareto scatter adds its own marker-size
  // multiselect earlier in the DOM, so a bare `.multiselect` would match that instead
  const select = document.querySelector<HTMLElement>(`.controls .multiselect`)
  if (!select) throw new Error(`No model multiselect found`)
  return select
}

const mount_page = async (search = ``): Promise<void> => {
  await mount_with_url(DiatomicsPage, `http://localhost/tasks/diatomics${search}`, {
    props: { data: page_data },
  })
}

const selected_options = (): HTMLElement | null =>
  model_select().querySelector(`ul[aria-label="selected options"]`)

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
  it(`uses default selected curve models without query params`, async () => {
    await mount_page()

    const selected_text = selected_options()?.textContent
    expect(selected_text).toContain(`PBE (DFT)`)
    expect(selected_text).toContain(`r2SCAN (DFT)`)
    expect(selected_text).toContain(`Model A`)
    expect(selected_text).toContain(`Model B`)
  })

  it(`restores selected curve models from the models query param`, async () => {
    await mount_page(`?models=model-b`)

    const selected = selected_options()
    expect(selected?.textContent).toContain(`Model B`)
    expect(selected?.textContent).not.toContain(`Model A`)
    expect(selected?.textContent).not.toContain(`PBE (DFT)`)
    expect(selected?.textContent).not.toContain(`r2SCAN (DFT)`)
    const selected_item = selected?.querySelector<HTMLLIElement>(`li`)
    expect(selected_item?.style.background).toContain(`#68d391`)
    expect(selected_item?.style.color).not.toBe(``)
  })

  it(`syncs selected curve models back to the models query param`, async () => {
    await mount_page(`?models=model-b`)

    await select_model_option(`Model A`)

    expect(new URL(location.href).searchParams.get(`models`)).toBe(`model-a,model-b`)
  })

  it(`restores selected element subset from the elements query param`, async () => {
    await mount_page(`?elements=halogen`)

    expect(button_for(`All`).getAttribute(`aria-pressed`)).toBe(`false`)
    expect(button_for(`Halogens`).getAttribute(`aria-pressed`)).toBe(`true`)
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
