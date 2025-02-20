import { mount, tick } from 'svelte'
import { beforeEach, describe, expect, test, vi } from 'vitest'
import { doc_query } from '.'
import Page from '../src/routes/tasks/diatomics/+page.svelte'

// Mock fetch function
const mockFetch = vi.fn()
vi.stubGlobal(`fetch`, mockFetch)

describe(`Diatomic Page`, () => {
  // Add container with dimensions to body before each test
  const container_style = `width: 800px; height: 600px;`
  beforeEach(() => {
    const container = document.createElement(`div`)
    container.setAttribute(`style`, container_style)
    document.body.appendChild(container)
    // Reset fetch mock
    mockFetch.mockReset()
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
    mockFetch.mockResolvedValueOnce({
      ok: true,
      json: () => Promise.resolve(sample_data),
    })

    mount(Page, { target: document.body })

    const title = doc_query(`h1`)
    expect(title.textContent).toContain(`Diatomic Energy Curves`)

    const controls = doc_query(`.controls`)
    expect(controls).toBeTruthy()
  })

  test(`loads and caches model data`, async () => {
    mockFetch.mockResolvedValue({
      ok: true,
      json: () => Promise.resolve(sample_data),
    })

    mount(Page, { target: document.body })

    // First load should make a fetch request for each model
    const n_initial_calls = mockFetch.mock.calls.length
    expect(n_initial_calls).toBe(0)

    // Trigger another load of the same model
    const button = doc_query(`button`) // click button to load another model
    button.click()
    await tick()

    // each model load first tries loading data from pred_file, then tries pred_file_url
    expect(mockFetch.mock.calls.length).toBe(n_initial_calls + 2)

    // Trigger another load of the same model
    button.click()
    await tick()

    // model should now be cached
    expect(mockFetch.mock.calls.length).toBe(n_initial_calls + 3)
  })

  test(`handles model toggle buttons`, async () => {
    mockFetch.mockResolvedValue({
      ok: true,
      json: () => Promise.resolve(sample_data),
    })

    mount(Page, { target: document.body })
    await tick()

    const buttons = document.querySelectorAll(`button`)
    expect(buttons.length).toBeGreaterThan(2)

    buttons[0].click()
    // Wait for data to load
    await tick()

    // Should show diatomic curves
    const plots = document.querySelectorAll(`.plot`)
    expect(plots.length).toBe(118)
  })

  test(`displays pretty model labels on buttons`, async () => {
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

  test(`applies correct button styling`, async () => {
    mount(Page, { target: document.body })

    const button = doc_query(`button`)
    const style = getComputedStyle(button)

    // Check button styling
    expect(style.borderWidth).toBe(`2px`)
    expect(style.background).toBe(`transparent`)
    expect(style.fontWeight).toBe(`500`)
  })
})
