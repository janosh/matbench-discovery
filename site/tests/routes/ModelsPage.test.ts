import { MODELS } from '$lib/models.svelte'
import { default as ModelsPage } from '$routes/models/+page.svelte'
import { mount, tick } from 'svelte'
import { describe, expect, it } from 'vitest'
import { doc_query } from '../index'

describe(`Models Page`, () => {
  it(`renders model sorting controls`, () => {
    mount(ModelsPage, { target: document.body })

    const toggle = document.querySelector(`input[type="checkbox"]`)
    expect(toggle?.parentElement?.textContent).toMatch(/show non-compliant models/i)

    expect(document.querySelector(`input[type="number"]`)).toBeDefined()
    expect(document.querySelectorAll(`input[type="radio"]`)).toHaveLength(2)
  })

  it(`renders metric sorting buttons`, () => {
    mount(ModelsPage, { target: document.body })

    const button_texts = Array.from(document.querySelectorAll(`ul button`)).map(
      (btn) => btn.textContent?.trim(),
    )
    expect(button_texts).toContain(`Model Name`)
    expect(button_texts).toContain(`F1`)
    expect(button_texts).toContain(`DAF`)
    expect(button_texts).toContain(`R2`)
  })

  it(`renders model cards`, () => {
    mount(ModelsPage, { target: document.body })

    const model_cards = document.querySelectorAll(`ol > li`)
    expect(model_cards.length).toBeGreaterThan(0)

    // Test first model card structure
    const first_card = model_cards[0]
    expect(first_card.querySelector(`h2 a`)).toBeDefined() // model name link
    expect(first_card.querySelector(`nav`)).toBeDefined() // links nav
    expect(first_card.querySelector(`.metrics`)).toBeDefined() // metrics section
  })

  it(`sorts models by selected metric`, async () => {
    mount(ModelsPage, { target: document.body })

    const initial_models = Array.from(document.querySelectorAll(`ol > li h2 a`)).map(
      (a) => a.textContent,
    )

    const daf_btn = document.querySelector(`button#DAF`) as HTMLButtonElement
    expect(daf_btn).toBeDefined()
    daf_btn?.click()
    await tick()

    const sorted_models = Array.from(document.querySelectorAll(`ol > li h2 a`)).map(
      (a) => a.textContent,
    )
    // Order should be different
    expect(sorted_models).not.toEqual(initial_models)
  })

  it(`toggles model details`, async () => {
    mount(ModelsPage, { target: document.body })

    const first_card = document.querySelector(`ol > li`)
    const details_btn = first_card?.querySelector(`h2 button`) as HTMLButtonElement
    expect(details_btn).toBeDefined()

    // Initially no authors section visible
    const initial_sections = first_card?.querySelectorAll(`section`)
    const has_authors = Array.from(initial_sections ?? []).some(
      (section) => section.querySelector(`h3`)?.textContent === `Authors`,
    )
    expect(has_authors).toBe(false)

    // Click to show details
    details_btn?.click()
    await tick()

    // Should now show authors and package versions
    const sections = first_card?.querySelectorAll(`section`)
    expect(sections?.length).toBeGreaterThan(0)
    expect(Array.from(sections ?? []).some((s) => s.textContent?.includes(`Authors`)))
      .toBe(
        true,
      )
    expect(
      Array.from(sections ?? []).some((s) => s.textContent?.includes(`Package versions`)),
    ).toBe(true)
  })

  it(`binds show_details state between page and model cards`, async () => {
    mount(ModelsPage, { target: document.body })

    const [first_card, second_card] = Array.from(document.querySelectorAll(`ol > li`))
    expect(first_card).toBeDefined()
    expect(second_card).toBeDefined()

    const first_details_btn = first_card.querySelector(`h2 button`) as HTMLButtonElement

    // Initially no details visible on either card
    // Check if the details sections are present by looking for non-metrics h3
    expect(first_card.querySelectorAll(`section:not(.metrics) h3`).length).toBe(0)
    expect(second_card.querySelectorAll(`section:not(.metrics) h3`).length).toBe(0)

    first_details_btn.click()
    await tick()

    // Now both cards should show details (shared state)
    // After clicking, details sections with h3 elements should be visible
    expect(first_card.querySelectorAll(`section:not(.metrics) h3`).length)
      .toBeGreaterThan(1)
    expect(second_card.querySelectorAll(`section:not(.metrics) h3`).length)
      .toBeGreaterThan(1)
  })

  it(`renders model limiting controls correctly`, () => {
    mount(ModelsPage, { target: document.body })

    const n_best_input = doc_query<HTMLInputElement>(`input[type="number"]`)

    expect(n_best_input).toBeDefined()
    expect(n_best_input.type).toBe(`number`)

    const initial_count = document.querySelectorAll(`ol > li`).length
    expect(initial_count).toBeGreaterThan(7)

    expect(n_best_input.min).toBe(`2`) // min_models = 2
    expect(n_best_input.max).toBe(String(initial_count))
    expect(Number(n_best_input.value)).toBe(initial_count)
    expect(document.querySelectorAll(`ol > li`).length).toBe(initial_count)

    const label_text = n_best_input.parentElement?.textContent
    expect(label_text).toMatch(/sort/i)
    expect(label_text).toMatch(/best models/i)
  })

  describe(`model limiting behavior`, () => {
    const min_models = 2

    it.each([
      { input: 3, expected: 3 },
      { input: 5, expected: 5 },
      { input: 10, expected: 10 },
      { input: MODELS.length - 1, expected: MODELS.length - 1 },
      { input: MODELS.length, expected: MODELS.length },
      { input: MODELS.length + 1, expected: MODELS.length }, // capped at max
      { input: MODELS.length + 100, expected: MODELS.length }, // far above max
      { input: 999, expected: MODELS.length }, // very large
      { input: 0, expected: min_models }, // below min
      { input: 1, expected: min_models }, // below min
      { input: -5, expected: min_models }, // negative
      { input: -999, expected: min_models }, // very negative
    ])(
      `displays $expected models when initial_show_n_best=$input`,
      ({ input, expected }) => {
        mount(ModelsPage, {
          target: document.body,
          props: { data: { initial_show_n_best: input } },
        })

        const model_cards = document.querySelectorAll(`ol > li`)
        const n_best_input = document.querySelector(
          `input[type="number"]`,
        ) as HTMLInputElement
        const effective = Math.min(Math.max(min_models, input), MODELS.length)

        expect(model_cards.length).toBe(expected)
        expect(Number(n_best_input.value)).toBe(effective)
        expect(n_best_input.min).toBe(String(min_models))
        expect(n_best_input.max).toBe(String(MODELS.length))
      },
    )

    it(`renders valid card structure with unique model names`, () => {
      mount(ModelsPage, {
        target: document.body,
        props: { data: { initial_show_n_best: 10 } },
      })

      const cards = document.querySelectorAll(`ol > li`)
      const names = Array.from(cards).map((card) => {
        expect(card.querySelector(`h2 a`)).toBeDefined()
        expect(card.querySelector(`.metrics`)).toBeDefined()
        const name = card.querySelector(`h2 a`)?.textContent?.trim()
        expect(name).toBeTruthy()
        return name
      })

      expect(new Set(names).size).toBe(names.length)
      expect(names.length).toBe(10)
    })
  })

  it(`maintains limit and structure when sorting changes`, async () => {
    const limit = 5
    mount(ModelsPage, {
      target: document.body,
      props: { data: { initial_show_n_best: limit } },
    })

    const n_best_input = document.querySelector(
      `input[type="number"]`,
    ) as HTMLInputElement
    expect(Number(n_best_input.value)).toBe(limit)

    const daf_btn = document.querySelector(`button#DAF`) as HTMLButtonElement
    daf_btn?.click()
    await tick()

    const cards = document.querySelectorAll(`ol > li`)
    expect(cards.length).toBeGreaterThan(0)
    expect(cards.length).toBeLessThanOrEqual(MODELS.length)
    cards.forEach((card) => {
      expect(card.querySelector(`h2 a`)).toBeDefined()
      expect(card.querySelector(`.metrics`)).toBeDefined()
    })
    expect(Number(n_best_input.value)).toBe(limit)
  })

  it(`shows correct subset when limiting models`, () => {
    mount(ModelsPage, {
      target: document.body,
      props: { data: { initial_show_n_best: 3 } },
    })

    const limited_names = Array.from(document.querySelectorAll(`ol > li h2 a`)).map((a) =>
      a.textContent
    )

    document.body.innerHTML = ``
    mount(ModelsPage, { target: document.body })

    const all_names = Array.from(document.querySelectorAll(`ol > li h2 a`)).map((a) =>
      a.textContent
    )
    expect(all_names.slice(0, 3)).toEqual(limited_names)
  })

  it(`defaults to showing all models with compliance filter checked`, () => {
    mount(ModelsPage, { target: document.body })

    expect(document.querySelectorAll(`ol > li`).length).toBe(MODELS.length)
    const checkbox = document.querySelector(`input[type="checkbox"]`) as HTMLInputElement
    expect(checkbox.checked).toBe(true)
  })

  it(`renders color legend`, () => {
    mount(ModelsPage, { target: document.body })

    const legend = document.querySelector(`legend`)
    expect(legend?.textContent).toContain(`best`)
    expect(legend?.textContent).toContain(`worst`)

    expect(legend?.querySelector(`svg`)).toBeDefined()
    expect(legend?.querySelector(`.matterviz-color-bar`)).toBeDefined()

    const model_cards_h2 = Array.from(
      document.querySelectorAll<HTMLElement>(`ol > li h2`),
    )
    expect(model_cards_h2.length).toBeGreaterThan(0)

    // applies background color to model card titles based on active metric value
    // currently only testing that the background color is not transparent
    for (const h2_element of model_cards_h2) {
      const computed_style = globalThis.getComputedStyle(h2_element)
      expect(computed_style.backgroundColor).not.toBe(`rgba(0, 0, 0, 0)`)
    }
  })
})
