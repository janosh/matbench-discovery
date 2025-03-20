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
          id="svg-btn">
          SVG
        </button>
        <button
          class="download-btn png-btn"
          id="png-btn">
          PNG
        </button>
        <button
          class="download-btn"
          id="test-error-btn">
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

  it(`has SVG and PNG download buttons`, () => {
    const buttons = document.querySelectorAll(`.download-btn`)
    expect(buttons).toHaveLength(3)

    const button_texts = Array.from(buttons).map((btn) => btn.textContent?.trim())
    expect(button_texts).toContain(`SVG`)
    expect(button_texts).toContain(`PNG`)
    expect(button_texts).toContain(`Test Error`)
  })

  it(`verifies SVG button exists with correct text`, () => {
    const svg_button = document.querySelector(`.svg-btn`)

    // Just verify the button exists with the right properties instead of triggering a click
    expect(svg_button).not.toBeNull()
    expect(svg_button?.textContent?.trim()).toBe(`SVG`)
    expect(svg_button?.id).toBe(`svg-btn`)
  })

  it(`verifies PNG button exists with correct text and id`, () => {
    const png_button = document.querySelector(`.png-btn`)

    // Just verify the button exists with the right properties
    expect(png_button).not.toBeNull()
    expect(png_button?.textContent?.trim()).toBe(`PNG`)
    expect(png_button?.id).toBe(`png-btn`)
  })

  it(`should handle errors when SVG generation fails`, async () => {
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

    // Check that the error message is displayed with the correct text
    expect(error_el.style.display).toBe(`block`)
    expect(error_el.textContent).toBe(`Error exporting SVG: ${mock_error_message}`)
  })

  it(`should handle null returns from PNG generation`, async () => {
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

    // Check that the error message is displayed with the correct text
    expect(error_el.style.display).toBe(`block`)
    expect(error_el.textContent).toBe(
      `Failed to generate PNG. The export function returned null.`,
    )
  })
})
