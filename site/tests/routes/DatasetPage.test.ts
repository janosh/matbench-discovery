import { arr_to_str, DATASETS } from '$lib'
import Page from '$routes/data/[slug]/+page.svelte'
import { beforeEach, describe, expect, it } from 'vitest'
import { doc_query, mount } from '../index'

describe(`Dataset Detail Page`, () => {
  beforeEach(() => {
    document.body.innerHTML = ``
  })

  // MP 2022 has all optional fields set, NOMAD is a minimal dataset entry
  it.each([`MP 2022`, `NOMAD`])(`renders %s dataset correctly`, (dataset_key) => {
    const dataset = DATASETS[dataset_key]
    if (!dataset) throw new Error(`Dataset ${dataset_key} not found in DATASETS`)

    mount(Page, { target: document.body, props: { data: { dataset } } })

    // Check title is displayed
    expect(document.querySelector(`h1`)?.textContent).toBe(dataset.name)

    // Check key metadata fields are displayed
    const meta_info = doc_query(`.meta-info`)
    expect(meta_info.textContent).toContain(`structures`)
    expect(meta_info.textContent).toContain(dataset.open ? `Open` : `Closed`)
    expect(meta_info.textContent).toContain(dataset.license)

    // Check links and description have content
    expect(doc_query(`.links`).querySelectorAll(`a`).length).toBeGreaterThan(0)
    expect(doc_query(`.description`).textContent).toMatch(/\S/)
  })

  // values used to be re-split on `:` after joining, truncating anything past a colon
  it(`lists method params with title-cased keys and colon-safe values`, () => {
    const dataset = DATASETS.WBM
    if (!dataset) throw new Error(`WBM not found in DATASETS`)
    const params = {
      code: `VASP`,
      cutoff_energy: `520 eV: hard`,
      pseudopotentials: [`PBE`],
    }
    mount(Page, {
      target: document.body,
      props: { data: { dataset: { ...dataset, params } } },
    })

    const items = [...doc_query(`.method-info`).querySelectorAll(`li`)].map((item) =>
      item.textContent?.replaceAll(/\s+/g, ` `).trim(),
    )
    // scalar values go through arr_to_str (JSON-quoted), arrays are comma-joined
    expect(items).toEqual([
      `Method: ${arr_to_str(dataset.method)}`,
      `Code: "VASP"`,
      `Cutoff Energy: "520 eV: hard"`,
      `Pseudopotentials: PBE`,
    ])
  })
})
