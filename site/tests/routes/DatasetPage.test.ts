import { DATASETS } from '$lib'
import Page from '$routes/data/[slug]/+page.svelte'
import { mount } from 'svelte'
import { beforeEach, describe, expect, it, vi } from 'vitest'

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
      props: { data: { dataset: mp_dataset, key: mp2022_key } },
    })

    // Check title is displayed
    expect(document.querySelector(`h1`)?.textContent).toBe(mp_dataset.title)

    // Check metadata section
    const meta_info = document.querySelector(`.meta-info`)
    expect(meta_info).not.toBeNull()

    // Check key metadata fields are displayed
    expect(meta_info?.textContent).toContain(`structures`)
    expect(meta_info?.textContent).toContain(mp_dataset.open ? `Open` : `Closed`)
    expect(meta_info?.textContent).toContain(mp_dataset.license)

    // Check links and description
    expect(document.querySelector(`.links`)).not.toBeNull()
    expect(document.querySelector(`.description`)).not.toBeNull()
  })

  it(`renders a minimal dataset correctly`, () => {
    // Skip test if no minimal dataset found
    if (!minimal_dataset) {
      console.warn(`Skipping test, couldn't find NOMAD dataset`)
      return
    }

    mount(Page, {
      target: document.body,
      props: { data: { dataset: minimal_dataset, key: minimal_key } },
    })

    // Check title and basic content
    expect(document.querySelector(`h1`)?.textContent).toBe(minimal_dataset.title)
    expect(document.querySelector(`.meta-info`)).not.toBeNull()
    expect(document.querySelector(`.description`)).not.toBeNull()

    // Check some specific fields
    const meta_info = document.querySelector(`.meta-info`)
    expect(meta_info?.textContent).toContain(`structures`)
    expect(meta_info?.textContent).toContain(minimal_dataset.open ? `Open` : `Closed`)
    expect(meta_info?.textContent).toContain(minimal_dataset.license)
  })
})
