import Page from '$routes/tasks/diatomics/+page.svelte'
import { mount } from 'svelte'
import { beforeEach, describe, expect, test, vi } from 'vitest'
import { doc_query } from '..'

// Mock fetch function
const mock_fetch = vi.fn()
vi.stubGlobal(`fetch`, mock_fetch)

describe(`Diatomic Page`, () => {
  // Add container with dimensions to body before each test
  const container_style = `width: 800px; height: 600px;`
  beforeEach(() => {
    const container = document.createElement(`div`)
    container.setAttribute(`style`, container_style)
    document.body.appendChild(container)
    // Reset fetch mock
    mock_fetch.mockReset()
  })

  const sample_data = {
    'homo-nuclear': {
      'H₂': {
        energies: [-5, -7, -3, -1, 0],
      },
    },
    distances: [0.5, 1.0, 1.5, 2.0, 2.5],
  }

  test(`renders page with initial state`, async () => {
    mock_fetch.mockResolvedValueOnce({
      ok: true,
      json: () => Promise.resolve(sample_data),
    })

    mount(Page, { target: document.body })

    const title = doc_query(`h1`)
    expect(title.textContent).toContain(`Diatomics`)

    const controls = doc_query(`.controls`)
    expect(controls).toBeTruthy()
  })

  test.skip(`loads and caches model data`, async () => {
    mock_fetch.mockResolvedValue({
      ok: true,
      json: () => Promise.resolve(sample_data),
    })

    mount(Page, { target: document.body })

    // First load should make a fetch request for each model
    const n_initial_calls = mock_fetch.mock.calls.length
    expect(n_initial_calls).toBe(0)

    // Trigger another load of the same model
    const button = doc_query(`button`) // click button to load another model
    button.click()

    // each model load first tries loading data from pred_file, then tries pred_file_url
    expect(mock_fetch.mock.calls.length).toBe(n_initial_calls + 2)

    button.click() // Trigger another load of the same model

    // model should now be cached
    expect(mock_fetch.mock.calls.length).toBe(n_initial_calls + 3)
  })

  test.skip(`displays pretty model labels on buttons`, async () => {
    mount(Page, { target: document.body })

    const buttons = document.querySelectorAll(`button`)
    const button_texts = Array.from(buttons).map((btn) => btn.textContent?.trim())

    // Check that buttons use pretty labels
    expect(button_texts).toContain(`MACE-MP-0`)
    // expect(button_texts).toContain(`MACE-MPA-0`)
  })

  test(`centers model toggle buttons`, async () => {
    mount(Page, { target: document.body })

    const controls = doc_query(`.controls`)
    const style = getComputedStyle(controls)
    expect(style.display).toBe(`flex`)
    expect(style[`flex-direction`]).toBe(`column`)
  })
})
