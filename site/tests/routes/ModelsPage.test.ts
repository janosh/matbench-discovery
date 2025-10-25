import { MODELS } from '$lib/models.svelte'
import { default as ModelsPage } from '$routes/models/+page.svelte'
import { mount, tick } from 'svelte'
import { afterEach, describe, expect, it } from 'vitest'

describe(`Models Page`, () => {
  afterEach(() => {
    document.body.innerHTML = ``
  })

  it(`renders model sorting controls`, () => {
    mount(ModelsPage, { target: document.body })

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
    mount(ModelsPage, { target: document.body })

    const buttons = document.querySelectorAll(`ul button`)
    const button_texts = Array.from(buttons).map((btn) => btn.textContent?.trim())

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
    expect(
      Array.from(sections ?? []).some((s) => s.textContent?.includes(`Authors`)),
    ).toBe(true)
    expect(
      Array.from(sections ?? []).some((s) => s.textContent?.includes(`Package versions`)),
    ).toBe(true)
  })

  it(`binds show_details state between page and model cards`, async () => {
    mount(ModelsPage, { target: document.body })

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
    mount(ModelsPage, { target: document.body })

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

  describe(`model limiting behavior`, () => {
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
      `limits displayed models to $expected when initial_show_n_best=$show_n_best`,
      ({ show_n_best, expected }) => {
        mount(ModelsPage, {
          target: document.body,
          props: { initial_show_n_best: show_n_best },
        })

        // Query the number input element
        const n_best_input = document.querySelector(
          `input[type="number"]`,
        ) as HTMLInputElement
        expect(n_best_input).toBeDefined()

        // Assert the rendered number of model cards
        const model_cards = document.querySelectorAll(`ol > li`)
        expect(model_cards.length).toBe(expected)

        // Verify the input shows the effective value (after Math.max)
        const effective_value = Math.max(min_models, show_n_best)
        expect(Number(n_best_input.value)).toBe(effective_value)
      },
    )
  })

  it(`renders color legend`, () => {
    mount(ModelsPage, { target: document.body })

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
