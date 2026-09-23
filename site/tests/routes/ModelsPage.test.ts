import { MODELS } from '$lib/models.svelte'
import { ALL_METRICS } from '$lib/labels'
import { sort_models } from '$lib/metrics'
import ModelsPage from '$routes/models/+page.svelte'
import { tick } from 'svelte'
import { describe, expect, it, vi } from 'vite-plus/test'
import { doc_query, mount, mount_with_url } from '../index'

// Outgoing cards become inert while their fade completes.
const model_cards = () => document.querySelectorAll(`ol.models > li:not([inert])`)

describe(`Models Page`, () => {
  it(`renders sorting controls and keyboard-accessible metric help`, async () => {
    mount(ModelsPage, { target: document.body })

    const n_best_input = doc_query<HTMLInputElement>(`input[type="number"]`)
    expect(n_best_input.value).toBe(String(MODELS.length))
    expect(n_best_input.min).toBe(`2`)
    expect(n_best_input.max).toBe(String(MODELS.length))
    expect(n_best_input.parentElement?.textContent).toMatch(/Sort\s+best models/)
    expect(
      [...document.querySelectorAll<HTMLInputElement>(`input[type="radio"]`)].map(
        (radio) => [radio.value, radio.checked],
      ),
    ).toStrictEqual([
      [`asc`, false],
      [`desc`, true],
    ])
    expect(
      [...document.querySelectorAll(`ul button`)].map((btn) => btn.textContent?.trim()),
    ).toStrictEqual([`Model Name`, `CPS`, `F1`, `RMSD`, `κSRME`, `CMDS`, `CDS`])
    const legend = doc_query(`legend`)
    expect(legend.textContent).toContain(`best`)
    expect(legend.textContent).toContain(`worst`)
    expect(doc_query(`.colorbar`, legend).textContent).toContain(
      `Card titles colored by CPS`,
    )
    for (const card of model_cards()) {
      expect(doc_query(`h2`, card).style.backgroundColor).not.toBe(``)
    }
    doc_query<HTMLButtonElement>(`button#CPS`).focus()
    await vi.waitFor(() =>
      expect(document.querySelector(`.popover`)?.textContent).toContain(
        `Combined Performance Score`,
      ),
    )
  })

  it(`renders every card with its canonical model link, including duplicate names`, () => {
    const original_name = MODELS[1].model_name
    MODELS[1].model_name = MODELS[0].model_name
    try {
      mount(ModelsPage, { target: document.body })
      const cards = model_cards()
      expect(cards).toHaveLength(MODELS.length)
      const expected_models = sort_models(MODELS, ALL_METRICS.CPS.key, `desc`)
      for (const [idx, card] of cards.entries()) {
        const { model_key, model_name } = expected_models[idx]
        const link = doc_query<HTMLAnchorElement>(`h2 a`, card)
        expect(link.getAttribute(`href`)).toBe(`/models/${model_key}`)
        expect(link.textContent?.trim()).toBe(model_name)
        expect(doc_query(`.metrics`, card).textContent).toContain(`CPS`)
      }
      expect(doc_query(`nav`, cards[0]).querySelectorAll(`a`).length).toBeGreaterThan(0)
      expect(cards[0].textContent).not.toContain(`Missing preds`)
      const activity = doc_query(`section[aria-labelledby="github-activity"]`)
      expect(doc_query(`a[href="#github-activity"]`).textContent).toContain(
        `GitHub activity`,
      )
      expect(doc_query(`#github-activity`, activity).textContent).toBe(`GitHub Activity`)
      expect(activity.querySelector(`div.scatter`)).not.toBeNull()
      expect(doc_query(`ol.models`).compareDocumentPosition(activity)).toBe(
        Node.DOCUMENT_POSITION_FOLLOWING,
      )
    } finally {
      MODELS[1].model_name = original_name
    }
  })

  it.each([
    [`F1`, `metrics.discovery.unique_prototypes.F1`, `desc`],
    [`diatomics_combined_score`, `metrics.diatomics.combined_score`, `desc`],
    [`rmsd`, `metrics.geo_opt.symprec=1e-2.rmsd`, `asc`],
    [`Model`, `model_name`, `asc`],
  ] as const)(
    `sorts by %s and preserves the model limit`,
    async (key, path, direction) => {
      await mount_with_url(ModelsPage, `http://localhost/models?n_best=5`)
      doc_query<HTMLButtonElement>(`button#${key}`).click()
      await tick()

      const expected_models = sort_models(MODELS, path, direction).slice(0, 5)
      const cards = model_cards()
      expect(
        [...cards].map((card) => doc_query(`h2 a`, card).getAttribute(`href`)),
      ).toStrictEqual(expected_models.map(({ model_key }) => `/models/${model_key}`))
      expect(doc_query(`ul li.active button`).id).toBe(key)
      expect(doc_query<HTMLInputElement>(`input[type="radio"]:checked`).value).toBe(
        direction,
      )
      expect(doc_query<HTMLInputElement>(`input[type="number"]`).value).toBe(`5`)
      for (const card of cards) {
        expect(card.querySelector(`.metrics li.active .metric-label`)?.innerHTML).toBe(
          Object.values(ALL_METRICS).find((metric) => metric.key === key)?.label,
        )
      }
      // Inspecting another metric must not activate the matching page-level sort button.
      doc_query(`.metrics .metric-label`).click()
      await tick()
      expect(doc_query(`ul li.active button`).id).toBe(key)
      const params = new URL(location.href).searchParams
      expect(params.get(`sort`)).toBe(key)
      expect(params.get(`dir`)).toBe(direction === `desc` ? null : direction)
      expect(params.get(`n_best`)).toBe(`5`)
    },
  )

  it.each([
    [``, MODELS.length],
    [`3.9`, 3],
    [`1e99`, MODELS.length],
    [`0`, MODELS.length],
    [`-3`, MODELS.length],
    [`abc`, MODELS.length],
  ])(`normalizes n_best=%s and invalid sorting parameters`, async (input, expected) => {
    await mount_with_url(
      ModelsPage,
      `http://localhost/models?sort=bogus&dir=sideways&n_best=${input}`,
    )

    expect(model_cards()).toHaveLength(expected)
    expect(doc_query<HTMLInputElement>(`input[type="number"]`).value).toBe(
      String(expected),
    )
    expect(doc_query(`ul li.active button`).id).toBe(`CPS`)
    expect(doc_query<HTMLInputElement>(`input[type="radio"]:checked`).value).toBe(`desc`)
  })

  it.each([
    [3, 3],
    [3.9, 3],
    [``, 2],
    [0, 2],
    [-5, 2],
    [MODELS.length + 1, MODELS.length],
  ] as const)(
    `restores URL state and applies a model limit of %s through the controls`,
    async (input, expected) => {
      await mount_with_url(ModelsPage, `http://localhost/models?sort=F1&dir=asc&n_best=5`)

      expect(doc_query(`ul li.active button`).id).toBe(`F1`)
      expect(doc_query<HTMLInputElement>(`input[type="radio"]:checked`).value).toBe(`asc`)
      expect(model_cards()).toHaveLength(5)
      const n_best_input = doc_query<HTMLInputElement>(`input[type="number"]`)
      expect(n_best_input.value).toBe(`5`)
      n_best_input.value = String(input)
      n_best_input.dispatchEvent(new Event(`input`, { bubbles: true }))
      doc_query<HTMLInputElement>(`input[type="radio"][value="desc"]`).click()
      await tick()

      expect(model_cards()).toHaveLength(expected)
      const params = new URL(location.href).searchParams
      expect(params.get(`sort`)).toBe(`F1`)
      expect(params.has(`dir`)).toBe(false)
      expect(params.get(`n_best`)).toBe(
        expected === MODELS.length ? null : String(expected),
      )
    },
  )
})
