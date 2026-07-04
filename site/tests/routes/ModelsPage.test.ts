import { MODELS } from '$lib/models.svelte'
import { ALL_METRICS, MD_METRICS } from '$lib/labels'
import { sort_models } from '$lib/metrics'
import { default as ModelsPage } from '$routes/models/+page.svelte'
import { mount, tick } from 'svelte'
import { describe, expect, it } from 'vitest'
import { doc_query, mount_with_url } from '../index'

describe(`Models Page`, () => {
  it(`renders model sorting controls`, () => {
    mount(ModelsPage, { target: document.body })

    const toggle = doc_query<HTMLInputElement>(`input[type="checkbox"]`)
    expect(toggle.checked).toBe(true)
    expect(toggle.parentElement?.textContent).toMatch(/show non-compliant models/i)

    const n_best_input = doc_query<HTMLInputElement>(`input[type="number"]`)
    expect(n_best_input.value).toBe(String(MODELS.length))
    expect(n_best_input.min).toBe(`2`) // min_models = 2
    expect(n_best_input.max).toBe(String(MODELS.length))
    const label_text = n_best_input.parentElement?.textContent
    expect(label_text).toMatch(/sort/i)
    expect(label_text).toMatch(/best models/i)
    expect(
      [...document.querySelectorAll<HTMLInputElement>(`input[type="radio"]`)].map(
        (radio) => [radio.value, radio.checked],
      ),
    ).toStrictEqual([
      // default sort direction is descending
      [`asc`, false],
      [`desc`, true],
    ])
  })

  it(`renders metric sorting buttons`, () => {
    mount(ModelsPage, { target: document.body })

    const button_texts = [...document.querySelectorAll(`ul button`)].map((btn) =>
      btn.textContent?.trim(),
    )
    // derive the full expected button order from the same label sources the page uses
    // (Model Name + 1-2 headline metrics per task), converting each label's HTML to
    // text the way the page's {@html label} render does, so removal or reorder of any
    // sort button is caught
    const html_to_text = (html: string): string => {
      const el = document.createElement(`div`)
      el.innerHTML = html
      return el.textContent?.trim() ?? ``
    }
    // mirrors metric_keys in site/src/routes/models/+page.svelte
    const headline_keys = [`CPS`, `F1`, `MAE`, `RMSD`, `κ_SRME`] as const
    const expected_button_texts = [
      `Model Name`,
      ...headline_keys.map((key) => html_to_text(ALL_METRICS[key].label)),
      html_to_text(MD_METRICS.md_combined_score.label),
      html_to_text(MD_METRICS.md_vdos_error.label),
    ]

    expect(button_texts).toStrictEqual(expected_button_texts)
  })

  it(`renders model cards`, () => {
    mount(ModelsPage, { target: document.body })

    const model_cards = document.querySelectorAll(`ol > li`)
    expect(model_cards).toHaveLength(MODELS.length)

    // Test first model card structure
    const first_card = model_cards[0]
    const first_link = doc_query<HTMLAnchorElement>(`h2 a`, first_card)
    expect(first_link.getAttribute(`href`)).toMatch(/^\/models\/[^/]+$/)
    expect(first_link.textContent?.trim()).toMatch(/\S/)
    expect(doc_query(`nav`, first_card).querySelectorAll(`a`).length).toBeGreaterThan(0)
    expect(doc_query(`.metrics`, first_card).textContent).toContain(`CPS`)
  })

  it(`sorts models by selected metric`, async () => {
    mount(ModelsPage, { target: document.body })

    const initial_models = [...document.querySelectorAll(`ol > li h2 a`)].map(
      (a) => a.textContent,
    )

    const f1_btn = doc_query<HTMLButtonElement>(`button#F1`)
    expect(f1_btn.textContent?.trim()).toBe(`F1`)
    f1_btn.click()
    await tick()

    const sorted_models = [...document.querySelectorAll(`ol > li h2 a`)].map(
      (a) => a.textContent,
    )
    // Order should be different
    expect(sorted_models).not.toStrictEqual(initial_models)
  })

  it(`toggles model details, shared across all cards via bound state`, async () => {
    mount(ModelsPage, {
      target: document.body,
      props: { data: { initial_show_n_best: 3 } },
    })

    const [first_card, second_card] = [...document.querySelectorAll(`ol > li`)]
    const details_btn = doc_query<HTMLButtonElement>(`h2 button`, first_card)
    expect(details_btn.getAttribute(`aria-label`)).toBe(
      `Show authors and package versions`,
    )
    expect(details_btn.getAttribute(`aria-expanded`)).toBe(`false`)

    // Initially no details sections (non-metrics h3s) visible on either card
    expect(first_card.querySelectorAll(`section:not(.metrics) h3`)).toHaveLength(0)
    expect(second_card.querySelectorAll(`section:not(.metrics) h3`)).toHaveLength(0)

    details_btn.click()
    await tick()
    expect(details_btn.getAttribute(`aria-expanded`)).toBe(`true`)

    const sections = [...first_card.querySelectorAll(`section`)]
    expect(sections.some((sec) => sec.textContent?.includes(`Authors`))).toBe(true)
    expect(sections.some((sec) => sec.textContent?.includes(`Package versions`))).toBe(
      true,
    )
    // show_details is bound page-wide, so the second card expands too
    expect(
      second_card.querySelectorAll(`section:not(.metrics) h3`).length,
    ).toBeGreaterThan(1)
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
        const n_best_input = doc_query<HTMLInputElement>(`input[type="number"]`)

        expect(model_cards).toHaveLength(expected)
        expect(Number(n_best_input.value)).toBe(expected)
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
      const expected_models = MODELS.toSorted(
        sort_models(ALL_METRICS.CPS.key, `desc`),
      ).slice(0, 10)
      const names = [...cards].map((card, idx) => {
        const link = doc_query<HTMLAnchorElement>(`h2 a`, card)
        expect(link.getAttribute(`href`)).toBe(
          `/models/${expected_models[idx].model_key}`,
        )
        expect(doc_query(`.metrics`, card).textContent).toContain(`CPS`)
        return link.textContent?.trim()
      })

      expect(names).toStrictEqual(expected_models.map((model) => model.model_name))
    })
  })

  it(`maintains limit and structure when sorting changes`, async () => {
    const limit = 5
    mount(ModelsPage, {
      target: document.body,
      props: { data: { initial_show_n_best: limit } },
    })

    const n_best_input = doc_query<HTMLInputElement>(`input[type="number"]`)
    expect(Number(n_best_input.value)).toBe(limit)

    const f1_btn = doc_query<HTMLButtonElement>(`button#F1`)
    f1_btn.click()
    await tick()

    const cards = document.querySelectorAll(`ol > li`)
    // Svelte transitions can briefly keep outgoing list items in the DOM after sorting.
    expect(cards.length).toBeGreaterThanOrEqual(limit)
    expect(cards.length).toBeLessThanOrEqual(MODELS.length)
    cards.forEach((card) => {
      expect(doc_query<HTMLAnchorElement>(`h2 a`, card).getAttribute(`href`)).toMatch(
        /^\/models\/[^/]+$/,
      )
      expect(doc_query(`.metrics`, card).textContent).toContain(`F1`)
    })
    expect(Number(n_best_input.value)).toBe(limit)
  })

  it(`shows correct subset when limiting models`, () => {
    mount(ModelsPage, {
      target: document.body,
      props: { data: { initial_show_n_best: 3 } },
    })

    const limited_names = [...document.querySelectorAll(`ol > li h2 a`)].map(
      (a) => a.textContent,
    )

    document.body.innerHTML = ``
    mount(ModelsPage, { target: document.body })

    const all_names = [...document.querySelectorAll(`ol > li h2 a`)].map(
      (a) => a.textContent,
    )
    expect(all_names.slice(0, 3)).toStrictEqual(limited_names)
  })

  it(`renders color legend`, () => {
    mount(ModelsPage, {
      target: document.body,
      props: { data: { initial_show_n_best: 3 } },
    })

    const legend = doc_query(`legend`)
    expect(legend.textContent).toContain(`best`)
    expect(legend.textContent).toContain(`worst`)

    expect(doc_query(`.colorbar`, legend).textContent).toContain(
      `Card titles colored by CPS`,
    )

    const model_cards_h2 = [...document.querySelectorAll<HTMLElement>(`ol > li h2`)]
    expect(model_cards_h2.length).toBeGreaterThan(0)

    // model_cards_h2 should receive inline backgroundColor styles from the page's
    // bg_color() computation; this does not validate transparency or exact colors.
    for (const h2_element of model_cards_h2) {
      expect(h2_element.style.backgroundColor).not.toBe(``)
    }
  })

  it(`restores sort key/dir, n_best and compliance filter from URL params`, async () => {
    await mount_with_url(
      ModelsPage,
      `http://localhost/models?sort=F1&dir=asc&n_best=5&non_compliant=0`,
    )

    expect(doc_query(`ul li.active button`).id).toBe(`F1`)
    expect(doc_query<HTMLInputElement>(`input[type="radio"]:checked`).value).toBe(`asc`)
    // card slicing from n_best is covered by the initial_show_n_best test; asserting
    // DOM counts here would race the out:fade transition on the re-render
    expect(doc_query<HTMLInputElement>(`input[type="number"]`).value).toBe(`5`)
    expect(doc_query<HTMLInputElement>(`input[type="checkbox"]`).checked).toBe(false)
  })

  it.each([
    `?sort=bogus&dir=sideways&n_best=1e99&non_compliant=maybe`,
    `?n_best=-3`,
    `?n_best=abc`,
  ])(`falls back to defaults for invalid URL params: %s`, async (query) => {
    await mount_with_url(ModelsPage, `http://localhost/models${query}`)

    expect(doc_query(`ul li.active button`).id).toBe(`CPS`)
    expect(doc_query<HTMLInputElement>(`input[type="radio"]:checked`).value).toBe(`desc`)
    expect(doc_query<HTMLInputElement>(`input[type="checkbox"]`).checked).toBe(true)
  })
})
