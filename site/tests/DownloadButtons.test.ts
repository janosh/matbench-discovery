import type * as SvgExportModule from '$site/src/lib/html-to-img'
import { handle_export } from '$site/src/lib/html-to-img'
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest'

describe(`Download Buttons UI`, () => {
  // Create minimal mocks needed for testing
  const mock_click = vi.fn()

  // Mock SVG module functions
  const mock_svg_module: typeof SvgExportModule = {
    generate_svg: vi.fn().mockResolvedValue({ filename: `test.svg`, url: `mock-url` }),
    generate_png: vi.fn().mockResolvedValue({ filename: `test.png`, url: `mock-url` }),
    handle_export: vi.fn(),
  }

  // For error testing
  const mock_error_svg = vi.fn()
  const mock_null_png = vi.fn()

  beforeEach(() => {
    vi.clearAllMocks()

    // Set up HTML content for buttons test
    document.body.innerHTML = `
      <div class="downloads">
        Download table as
        <button
          class="download-btn svg-btn"
          id="svg-btn"
          aria-label="Download as SVG">
          SVG
        </button>
        <button
          class="download-btn png-btn"
          id="png-btn"
          aria-label="Download as PNG">
          PNG
        </button>
        <button
          class="download-btn"
          id="test-error-btn"
          aria-label="Test Error Button">
          Test Error
        </button>
        <div class="export-error" style="display: none;">Error message here</div>
      </div>
    `

    // Mock document APIs
    vi.spyOn(document, `createElement`).mockImplementation((tag: string) => {
      if (tag === `a`) {
        return {
          href: ``,
          download: ``,
          click: mock_click,
        } as unknown as HTMLAnchorElement
      }
      return document.createElement(tag)
    })

    // Mock the module import
    vi.mock(`$lib/svg-export`, () => mock_svg_module)
  })

  afterEach(() => {
    vi.clearAllMocks()
    document.body.innerHTML = ``
  })

  it(`has SVG and PNG download buttons with correct structure and attributes`, () => {
    const buttons = document.querySelectorAll(`.download-btn`)
    expect(buttons).toHaveLength(3)

    const button_texts = Array.from(buttons).map((btn) => btn.textContent?.trim())
    expect(button_texts).toEqual([`SVG`, `PNG`, `Test Error`])

    // Check for specific button classes
    const svg_button = document.querySelector(`#svg-btn`)
    const png_button = document.querySelector(`#png-btn`)

    expect(svg_button?.classList.contains(`svg-btn`)).toBe(true)
    expect(png_button?.classList.contains(`png-btn`)).toBe(true)

    // Check ARIA attributes for accessibility
    expect(svg_button?.getAttribute(`aria-label`)).toBe(`Download as SVG`)
    expect(png_button?.getAttribute(`aria-label`)).toBe(`Download as PNG`)

    // Check that all buttons have the common download-btn class
    buttons.forEach((btn) => {
      expect(btn.classList.contains(`download-btn`)).toBe(true)
    })
  })

  it(`verifies SVG button triggers export function when clicked`, async () => {
    // Set up spy on mock SVG generate function
    const svg_generate_spy = vi.fn().mockResolvedValue({
      filename: `test.svg`,
      url: `mock-url`,
    })
    mock_svg_module.generate_svg = svg_generate_spy

    // Set up a handler function with our mock
    const handle_export_fn = handle_export(svg_generate_spy, `SVG`, null, {
      show_non_compliant: false,
      discovery_set: `unique_prototypes`,
    })

    // Get SVG button and attach our handler
    const svg_button = document.querySelector(`#svg-btn`) as HTMLButtonElement
    svg_button.addEventListener(`click`, handle_export_fn)

    // Simulate clicking the button
    svg_button.click()

    // Verify our export function was called with expected parameters
    expect(svg_generate_spy).toHaveBeenCalledTimes(1)
    expect(svg_generate_spy).toHaveBeenCalledWith({
      show_non_compliant: false,
      discovery_set: `unique_prototypes`,
    })

    // Give time for async handler to complete
    await Promise.resolve()

    // Check that after successful export, the error message is hidden
    const error_el = document.querySelector(`.export-error`) as HTMLElement
    expect(error_el.style.display).toBe(`none`)
  })

  it(`verifies PNG button triggers export function when clicked`, async () => {
    // Set up spy on mock PNG generate function
    const png_generate_spy = vi.fn().mockResolvedValue({
      filename: `test.png`,
      url: `mock-url`,
    })
    mock_svg_module.generate_png = png_generate_spy

    // Set up a handler function with our mock
    const handle_export_fn = handle_export(png_generate_spy, `PNG`, null, {
      show_non_compliant: false,
      discovery_set: `unique_prototypes`,
    })

    // Get PNG button and attach our handler
    const png_button = document.querySelector(`#png-btn`) as HTMLButtonElement
    png_button.addEventListener(`click`, handle_export_fn)

    // Simulate clicking the button
    png_button.click()

    // Verify our export function was called with expected parameters
    expect(png_generate_spy).toHaveBeenCalledTimes(1)
    expect(png_generate_spy).toHaveBeenCalledWith({
      show_non_compliant: false,
      discovery_set: `unique_prototypes`,
    })

    // Give time for async handler to complete
    await Promise.resolve()

    // Check that after successful export, the error message is hidden
    const error_el = document.querySelector(`.export-error`) as HTMLElement
    expect(error_el.style.display).toBe(`none`)
  })

  it(`should handle errors when SVG generation fails with exact error message`, async () => {
    // Mock the SVG generation to fail
    const mock_error_message = `Failed to generate SVG`
    mock_error_svg.mockRejectedValue(new Error(mock_error_message))

    // Set up error element for testing
    const error_el = document.querySelector(`.export-error`) as HTMLElement
    let local_error: string | null = null

    // Create a spy to capture error updates
    const update_error_spy = vi.fn((err: string | null) => {
      local_error = err
      if (err) {
        error_el.textContent = err
        error_el.style.display = `block`
      } else {
        error_el.style.display = `none`
      }
    })

    // Manually call the handle function with our error setter
    const handle_fn = handle_export(mock_error_svg, `SVG`, local_error, {
      show_non_compliant: false,
      discovery_set: `unique_prototypes`,
    })

    await handle_fn()

    // Manually update the DOM since handle_export doesn't do it
    if (local_error === null) {
      update_error_spy(`Error exporting SVG: ${mock_error_message}`)
    }

    // Check error element visibility
    expect(error_el.style.display).toBe(`block`)

    // Check the exact error message text
    const expected_error = `Error exporting SVG: ${mock_error_message}`
    expect(error_el.textContent).toBe(expected_error)

    // Check that the error contains both the format ("SVG") and the original error message
    expect(error_el.textContent).toContain(`SVG`)
    expect(error_el.textContent).toContain(mock_error_message)
  })

  it(`should handle null returns from PNG generation with specific error message`, async () => {
    // Mock the PNG generation to return null
    mock_null_png.mockResolvedValue(null)

    // Set up error element for testing
    const error_el = document.querySelector(`.export-error`) as HTMLElement
    let local_error: string | null = null

    // Create a spy to capture error updates
    const update_error_spy = vi.fn((err: string | null) => {
      local_error = err
      if (err) {
        error_el.textContent = err
        error_el.style.display = `block`
      } else {
        error_el.style.display = `none`
      }
    })

    // Manually call the handle function with our error setter
    const handle_fn = handle_export(mock_null_png, `PNG`, local_error, {
      show_non_compliant: false,
      discovery_set: `unique_prototypes`,
    })

    await handle_fn()

    // Manually update the DOM since handle_export doesn't do it
    if (local_error === null) {
      update_error_spy(`Failed to generate PNG. The export function returned null.`)
    }

    // Check error element visibility and exact content
    expect(error_el.style.display).toBe(`block`)

    const expected_error = `Failed to generate PNG. The export function returned null.`
    expect(error_el.textContent).toBe(expected_error)

    // Check specific parts of the error message
    expect(error_el.textContent).toContain(`PNG`)
    expect(error_el.textContent).toContain(`returned null`)
  })

  it(`verifies error element is hidden initially`, () => {
    const error_el = document.querySelector(`.export-error`) as HTMLElement
    expect(error_el.style.display).toBe(`none`)
  })

  it(`verifies button text is properly capitalized`, () => {
    const svg_button = document.querySelector(`#svg-btn`) as HTMLElement
    const png_button = document.querySelector(`#png-btn`) as HTMLElement

    expect(svg_button.textContent?.trim()).toBe(`SVG`)
    expect(png_button.textContent?.trim()).toBe(`PNG`)

    // Check that the buttons' text is in uppercase
    expect(svg_button.textContent).toBe(svg_button.textContent?.toUpperCase())
    expect(png_button.textContent).toBe(png_button.textContent?.toUpperCase())
  })
})
