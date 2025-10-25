import { MODELS } from '$lib/models.svelte'
import { default as ModelsPage } from '$routes/models/+page.svelte'
import { mount, tick } from 'svelte'
import { beforeEach, describe, expect, it } from 'vitest'

describe(`Models Page`, () => {
  beforeEach(() => {
    mount(ModelsPage, { target: document.body })
  })

  it(`renders model sorting controls`, () => {
    const toggle = document.querySelector(`input[type="checkbox"]`)
    expect(toggle).toBeDefined()
    expect(toggle?.parentElement?.textContent).toMatch(/show non-compliant models/i)

    const n_best_input = document.querySelector(
      `input[type="number"]`,
    ) as HTMLInputElement
    expect(n_best_input).toBeDefined()

    const radio_buttons = document.querySelectorAll(`input[type="radio"]`)
    expect(radio_buttons).toHaveLength(2)
  })

  it(`renders metric sorting buttons`, () => {
    const buttons = document.querySelectorAll(`ul button`)
    const button_texts = Array.from(buttons).map((btn) => btn.textContent?.trim())

    expect(button_texts).toContain(`Model Name`)
    expect(button_texts).toContain(`F1`)
    expect(button_texts).toContain(`DAF`)
    expect(button_texts).toContain(`R2`)
  })

  it(`renders model cards`, () => {
    const model_cards = document.querySelectorAll(`ol > li`)
    expect(model_cards.length).toBeGreaterThan(0)

    // Test first model card structure
    const first_card = model_cards[0]
    expect(first_card.querySelector(`h2 a`)).toBeDefined() // model name link
    expect(first_card.querySelector(`nav`)).toBeDefined() // links nav
    expect(first_card.querySelector(`.metrics`)).toBeDefined() // metrics section
  })

  it(`sorts models by selected metric`, async () => {
    // Get initial order of models
    const initial_models = Array.from(document.querySelectorAll(`ol > li h2 a`)).map(
      (a) => a.textContent,
    )

    // Click DAF button to sort by DAF
    const daf_btn = document.querySelector(`button#DAF`) as HTMLButtonElement
    expect(daf_btn).toBeDefined()
    daf_btn?.click()
    await tick() // wait for re-render

    // Get new order of models
    const sorted_models = Array.from(document.querySelectorAll(`ol > li h2 a`)).map((a) =>
      a.textContent
    )

    // Order should be different
    expect(sorted_models).not.toEqual(initial_models)
  })

  it(`toggles model details`, async () => {
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
    expect(
      Array.from(sections ?? []).some((s) => s.textContent?.includes(`Authors`)),
    ).toBe(true)
    expect(
      Array.from(sections ?? []).some((s) => s.textContent?.includes(`Package versions`)),
    ).toBe(true)
  })

  it(`binds show_details state between page and model cards`, async () => {
    const model_cards = Array.from(
      document.querySelectorAll<HTMLElement>(`ol > li`),
    )
    expect(model_cards.length).toBeGreaterThan(1)

    const [first_card, second_card] = model_cards
    const first_details_btn = first_card.querySelector(`h2 button`) as HTMLButtonElement

    // Initially no details visible on either card
    // Check if the details sections are present by looking for non-metrics h3
    expect(first_card.querySelectorAll(`section:not(.metrics) h3`).length).toBe(0)
    expect(second_card.querySelectorAll(`section:not(.metrics) h3`).length).toBe(0)

    first_details_btn.click()
    await tick()
    // Now both cards should show details (shared state)
    // After clicking, details sections with h3 elements should be visible
    expect(
      first_card.querySelectorAll(`section:not(.metrics) h3`).length,
    ).toBeGreaterThan(1)
    expect(
      second_card.querySelectorAll(`section:not(.metrics) h3`).length,
    ).toBeGreaterThan(1)
  })

  it(`renders model limiting controls correctly`, () => {
    const n_best_input = document.querySelector(
      `input[type="number"]`,
    ) as HTMLInputElement

    expect(n_best_input).toBeDefined()
    expect(n_best_input.type).toBe(`number`)

    const initial_count = document.querySelectorAll(`ol > li`).length
    expect(initial_count).toBeGreaterThan(7)

    // Verify the input is properly constrained
    expect(n_best_input.min).toBe(`2`) // min_models = 2
    expect(n_best_input.max).toBe(String(initial_count))
    expect(Number(n_best_input.value)).toBe(initial_count)
    expect(document.querySelectorAll(`ol > li`).length).toBe(initial_count)

    const label_text = n_best_input.parentElement?.textContent
    expect(label_text).toMatch(/sort/i)
    expect(label_text).toMatch(/best models/i)
  })

  describe(`slice logic`, () => {
    const min_models = 2

    it.each([
      { show_n_best: 3, expected: 3 },
      { show_n_best: 5, expected: 5 },
      { show_n_best: 10, expected: 10 },
      { show_n_best: MODELS.length, expected: MODELS.length },
      { show_n_best: 0, expected: min_models },
      { show_n_best: 1, expected: min_models },
      { show_n_best: -5, expected: min_models },
    ])(
      `limits displayed models with show_n_best=$show_n_best`,
      ({ show_n_best, expected }) => {
        const limit = Math.max(min_models, show_n_best)
        const sliced = MODELS.slice(0, limit)
        expect(sliced.length).toBe(expected)
      },
    )

    it(`enforces minimum of ${min_models} models with invalid inputs`, () => {
      const invalid_inputs = [0, -1, -100, Number.NEGATIVE_INFINITY]
      for (const show_n_best of invalid_inputs) {
        expect(Math.max(min_models, show_n_best)).toBe(min_models)
      }
    })

    it(`respects exact number when above minimum`, () => {
      const valid_inputs = [3, 5, 10, 20, MODELS.length]
      for (const show_n_best of valid_inputs) {
        const sliced = MODELS.slice(0, Math.max(min_models, show_n_best))
        expect(sliced.length).toBe(Math.min(show_n_best, MODELS.length))
      }
    })
  })

  it(`renders color legend`, () => {
    const legend = document.querySelector(`legend`)
    expect(legend?.textContent).toContain(`best`)
    expect(legend?.textContent).toContain(`worst`)

    // Check that the ColorBar component rendered its SVG
    const color_bar_svg = legend?.querySelector(`svg`)
    expect(color_bar_svg).toBeDefined()

    const color_bar = legend?.querySelector(`.matterviz-color-bar`)
    expect(color_bar).toBeDefined()

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
