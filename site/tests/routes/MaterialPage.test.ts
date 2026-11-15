import kappa_analysis from '#figs/kappa-103-analysis.jsonl'
import { kappa_material } from '#lib/materials.js'
import { has_kappa_modes } from '#lib/parity/kappa-modes.js'
import type * as KappaModes from '#lib/parity/kappa-modes.js'
import type * as KappaParity from '#lib/parity/kappa-parity.js'
import { has_kappa_parity_model } from '#lib/parity/kappa-parity.js'
import MaterialPage from '#routes/materials/phonons/[material_id]/+page.svelte'
import { describe, expect, it, vi } from 'vite-plus/test'
import { heading_texts, header_name, mount, selected_chips } from '../index.js'

// asset loads never settle: the test checks what the page requests, not the plots
const { loads, pending } = vi.hoisted(() => ({
  loads: { parity: [] as string[], modes: [] as string[] },
  pending: () => new Promise<never>(() => {}),
}))
vi.mock(`#lib/parity/kappa-parity.js`, async (import_original) => ({
  ...(await import_original<typeof KappaParity>()),
  load_kappa_parity_base: pending,
  load_kappa_parity_model: (model_key: string) => {
    loads.parity.push(model_key)
    return pending()
  },
}))
vi.mock(`#lib/parity/kappa-modes.js`, async (import_original) => ({
  ...(await import_original<typeof KappaModes>()),
  load_kappa_modes: (model_key: string) => {
    loads.modes.push(model_key)
    return pending()
  },
}))

describe(`material page`, () => {
  it(`compares all models and loads the DOS and band structure of the right ones`, () => {
    const material = kappa_material(kappa_analysis, `mp-1184046`)
    const { best, worst, rows } = material
    if (!best || !worst) throw new Error(`fixture material has no scored models`)
    mount(MaterialPage, { target: document.body, props: { data: { material } } })

    expect(document.querySelectorAll(`tbody tr`)).toHaveLength(rows.length)
    expect(heading_texts()).not.toContain(`PBE reference structure`)
    // DOS defaults to the best and worst scored model that has a parity asset
    const dos_defaults = [best, worst]
      .filter(({ model_key }) => has_kappa_parity_model(model_key))
      .map(({ label }) => label)
    expect(selected_chips(`.dos-models`)).toEqual(dos_defaults)
    expect(loads.parity).toEqual([best, worst].map(({ model_key }) => model_key))
    // the band structure loads up front, for the best-ranked model with a modes asset
    const first_with_modes = rows.find(({ model_key }) => has_kappa_modes(model_key))
    expect(loads.modes).toEqual([first_with_modes?.model_key])
  })

  it.each([
    [`mp-1184046`, false],
    [`mp-24093`, true],
  ])(
    `%s: Issues column shown only when a model failed (%s)`,
    (material_id, has_issues) => {
      const material = kappa_material(kappa_analysis, material_id)
      mount(MaterialPage, { target: document.body, props: { data: { material } } })
      const headers = [...document.querySelectorAll(`thead th`)].map(header_name)
      expect(headers.includes(`Issues`)).toBe(has_issues)
      // the censored-marker note only when such markers exist
      expect(document.body.textContent?.includes(`Hollow red markers`)).toBe(has_issues)
      // prev/next navigation sits below the content, not above the title
      expect(document.querySelector(`section ~ .material-nav`)).not.toBeNull()
      const nav_hrefs = [...document.querySelectorAll(`.material-nav a`)].map((link) =>
        link.getAttribute(`href`),
      )
      expect(nav_hrefs).toContain(`/materials/phonons`)
      expect(nav_hrefs.every((href) => href?.startsWith(`/materials/phonons`))).toBe(true)
      expect(
        document.querySelector(`.material-nav ~ section, .material-nav ~ h1`),
      ).toBeNull()
      expect(document.body.textContent).not.toContain(`Materials Project`)
      const facts = document.querySelector(`.facts`)?.textContent ?? ``
      expect(facts).toContain(`#${material.median_srme_rank}/${material.n_ranked}`)
      expect(facts).toMatch(/Within 2× · 30% · 10% of DFT κ\s*\d+ · \d+ · \d+/)
    },
  )
})
