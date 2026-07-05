import { DATASETS } from '$lib'
import Page from '$routes/data/[slug]/+page.svelte'
import { beforeEach, describe, expect, it, vi } from 'vitest'
import { doc_query, mount } from '../index'

// Find two real datasets to test with
const mp2022_key = `MP 2022`
const minimal_key = `NOMAD`
// Get actual datasets
const mp_dataset = DATASETS[mp2022_key]
const minimal_dataset = DATASETS[minimal_key]

describe(`Dataset Detail Page`, () => {
  beforeEach(() => {
    document.body.innerHTML = ``
    vi.clearAllMocks()
  })

  it(`renders a dataset correctly with all optional fields`, () => {
    mount(Page, {
      target: document.body,
      props: { data: { dataset: mp_dataset } },
    })

    // Check title is displayed
    expect(document.querySelector(`h1`)?.textContent).toBe(mp_dataset.name)

    // Check key metadata fields are displayed
    const meta_info = doc_query(`.meta-info`)
    expect(meta_info.textContent).toContain(`structures`)
    expect(meta_info.textContent).toContain(mp_dataset.open ? `Open` : `Closed`)
    expect(meta_info.textContent).toContain(mp_dataset.license)

    // Check links and description have content
    expect(doc_query(`.links`).querySelectorAll(`a`).length).toBeGreaterThan(0)
    expect(doc_query(`.description`).textContent).toMatch(/\S/)
  })

  it(`renders a minimal dataset correctly`, () => {
    // Skip test if no minimal dataset found
    if (!minimal_dataset) {
      console.warn(`Skipping test, couldn't find NOMAD dataset`)
      return
    }

    mount(Page, {
      target: document.body,
      props: { data: { dataset: minimal_dataset } },
    })

    // Check title and basic content
    expect(document.querySelector(`h1`)?.textContent).toBe(minimal_dataset.name)
    expect(doc_query(`.description`).textContent).toMatch(/\S/)

    // Check some specific fields
    const meta_info = doc_query(`.meta-info`)
    expect(meta_info.textContent).toContain(`structures`)
    expect(meta_info.textContent).toContain(minimal_dataset.open ? `Open` : `Closed`)
    expect(meta_info.textContent).toContain(minimal_dataset.license)
  })
})
