import { ACTIVE_MODELS } from '$lib/models.svelte'
import ElementErrorsPtableHeatmap from '$routes/tasks/discovery/tmi/ElementErrorsPtableHeatmap.svelte'
import { per_element_each_errors as per_elem_each_errors } from '$lib/per-element-errors'
import { describe, expect, it } from 'vitest'
import { tick } from 'svelte'
import { format_num } from 'matterviz/labels'
import { doc_query, mount_with_url } from '../index'

const models_with_errors = ACTIVE_MODELS.filter(
  ({ model_key }) => model_key in per_elem_each_errors,
)

describe(`ElementErrorsPtableHeatmap`, () => {
  it.each([
    [`alchembert`, false],
    [`alchembert`, true],
    [`alchembert,chgnet-0.3.0`, false],
    [`alchembert,chgnet-0.3.0`, true],
  ] as const)(
    `distinguishes missing errors and units for %s, normalized=%s`,
    async (keys, normalized) => {
      await mount_with_url(
        ElementErrorsPtableHeatmap,
        `http://localhost/tasks/discovery/tmi?element_models=${keys}&element_normalized=${Number(normalized)}`,
      )
      const helium = doc_query(`[data-element-symbol="He"]`)
      helium.dispatchEvent(new MouseEvent(`mouseenter`))
      await tick()
      const inset = doc_query(`.model-errors`)
      expect(inset.textContent).toContain(`Helium:`)
      expect([...inset.querySelectorAll(`b`)].map((value) => value.textContent)).toEqual(
        keys.split(`,`).map(() => `n/a`),
      )
      expect(helium.querySelector(`.value`)).toBeNull()
      expect(helium.style.backgroundColor).toBe(`rgba(255, 255, 255, 0.3)`)
      expect(doc_query(`small`, inset).textContent).toBe(
        normalized ? `normalized` : `eV/atom`,
      )

      doc_query(`[data-element-symbol="H"]`).dispatchEvent(new MouseEvent(`mouseenter`))
      await tick()
      const std = normalized ? 0.1709639053 : 1
      expect([...inset.querySelectorAll(`b`)].map((value) => value.textContent)).toEqual(
        (keys.includes(`,`) ? [0.5663, 0.3445] : [0.5663]).map((value) =>
          format_num(value / std),
        ),
      )
    },
  )

  it.each([``, `?element_models=unknown&element_max=Infinity`])(
    `defaults to a model with per-element error data: %s`,
    async (query) => {
      await mount_with_url(
        ElementErrorsPtableHeatmap,
        `http://localhost/tasks/discovery/tmi${query}`,
      )

      const chips = document.querySelectorAll(`ul[aria-label="selected options"] li`)
      expect(chips).toHaveLength(1)
      expect(chips[0].textContent?.trim()).toBe(models_with_errors[0].model_name)
      expect(location.search).toBe(``)
    },
  )

  it(`restores ordered model selections and writes heatmap controls to the URL`, async () => {
    const models = models_with_errors.slice(-5).toReversed()
    const keys = models.map(({ model_key }) => model_key)
    await mount_with_url(
      ElementErrorsPtableHeatmap,
      `http://localhost/tasks/discovery/tmi?element_models=${keys[0]},unknown,${keys.join(`,`)}&element_normalized=0&element_manual_max=1&element_max=0.4`,
    )
    const chips = [...document.querySelectorAll(`ul[aria-label="selected options"] li`)]
    expect(chips.map((chip) => chip.textContent?.trim())).toEqual(
      models.slice(0, 4).map(({ model_name }) => model_name),
    )
    expect(new URLSearchParams(location.search).get(`element_models`)).toBe(
      keys.slice(0, 4).join(`,`),
    )
    const [manual, normalized] = document.querySelectorAll<HTMLInputElement>(
      `form input[type="checkbox"]`,
    )
    const maximum = doc_query<HTMLInputElement>(`input[type="range"]`)
    expect(manual.checked).toBe(true)
    expect(normalized.checked).toBe(false)
    expect(maximum.value).toBe(`0.4`)
    maximum.value = `0.5`
    maximum.dispatchEvent(new Event(`input`, { bubbles: true }))
    normalized.click()
    manual.click()
    doc_query<HTMLButtonElement>(`button`, chips[1]).click()
    await tick()
    const params = new URLSearchParams(location.search)
    expect(params.get(`element_max`)).toBe(`0.5`)
    expect(params.get(`element_models`)).toBe([keys[0], keys[2], keys[3]].join(`,`))
    expect(params.has(`element_manual_max`)).toBe(false)
    expect(params.has(`element_normalized`)).toBe(false)
  })
})
