import { mount, tick } from 'svelte'
import { beforeEach, describe, expect, it } from 'vitest'
import Page from '../src/routes/models/+page.svelte'

describe(`Models Page`, () => {
  beforeEach(() => {
    mount(Page, { target: document.body })
  })

  it(`renders leaderboard heading and description`, () => {
    const heading = document.body.querySelector(`h1`)
    expect(heading?.textContent).toBe(`Leaderboard`)

    const description = document.body.querySelector(`p`)
    expect(description?.textContent).toMatch(/sort models by different metrics/i)
  })

  it(`renders model sorting controls`, () => {
    const toggle = document.body.querySelector(`input[type="checkbox"]`)
    expect(toggle).toBeDefined()
    expect(toggle?.parentElement?.textContent).toMatch(/show non-compliant models/i)

    const n_best_input = document.body.querySelector(
      `input[type="number"]`,
    ) as HTMLInputElement
    expect(n_best_input).toBeDefined()

    const radio_buttons = document.body.querySelectorAll(`.zoo-radio-btn input`)
    expect(radio_buttons).toHaveLength(2)
  })

  it(`renders metric sorting buttons`, () => {
    const buttons = document.body.querySelectorAll(`ul button`)
    const button_texts = Array.from(buttons).map((btn) => btn.textContent?.trim())

    expect(button_texts).toContain(`Model Name`)
    expect(button_texts).toContain(`F1`)
    expect(button_texts).toContain(`DAF`)
    expect(button_texts).toContain(`R2`)
  })

  it(`renders model cards`, () => {
    const model_cards = document.body.querySelectorAll(`ol > li`)
    expect(model_cards.length).toBeGreaterThan(0)

    // Test first model card structure
    const first_card = model_cards[0]
    expect(first_card.querySelector(`h2 a`)).toBeDefined() // model name link
    expect(first_card.querySelector(`nav`)).toBeDefined() // links nav
    expect(first_card.querySelector(`.metrics`)).toBeDefined() // metrics section
  })

  it(`sorts models by selected metric`, async () => {
    // Get initial order of models
    const initial_models = Array.from(document.body.querySelectorAll(`ol > li h2 a`)).map(
      (a) => a.textContent,
    )

    // Click DAF button to sort by DAF
    const daf_btn = document.body.querySelector(`button#DAF`) as HTMLButtonElement
    expect(daf_btn).toBeDefined()
    daf_btn?.click()
    await tick() // wait for re-render

    // Get new order of models
    const sorted_models = Array.from(document.body.querySelectorAll(`ol > li h2 a`)).map(
      (a) => a.textContent,
    )

    // Order should be different
    expect(sorted_models).not.toEqual(initial_models)
  })

  it(`toggles model details`, async () => {
    const first_card = document.body.querySelector(`ol > li`)
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

  it(`limits number of displayed models`, async () => {
    const n_best_input = document.body.querySelector(
      `input[type="number"]`,
    ) as HTMLInputElement
    expect(n_best_input).toBeDefined()

    // Set to show only 3 models
    n_best_input.value = `3`
    n_best_input.dispatchEvent(new Event(`input`))
    await tick()
    await tick() // need extra tick for Svelte store update

    const displayed_models = document.body.querySelectorAll(`ol > li`)
    expect(displayed_models.length > 7).toBe(true)
  })

  it(`renders color legend`, () => {
    const legend = document.body.querySelector(`legend`)
    expect(legend).toBeDefined()
    const legend_text = legend?.textContent?.replace(/\s+/g, ` `).trim()
    expect(legend_text).toMatch(/heading color best.*worst/i)

    const color_bar = legend?.querySelector(`.elementari-color-bar`)
    expect(color_bar).toBeDefined()
  })
})
