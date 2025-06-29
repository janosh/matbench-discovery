import { handle_export } from '$lib/table-export'
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest'

type ExportState = {
  export_error: string | null
  show_non_compliant: boolean
  discovery_set: string
}

describe(`Download Buttons UI`, () => {
  const mock_click = vi.fn()

  const default_state: ExportState = {
    export_error: null,
    show_non_compliant: false,
    discovery_set: `unique_prototypes`,
  }

  beforeEach(() => {
    vi.clearAllMocks()
    document.body.innerHTML = `
      <div class="downloads">
        <button class="download-btn" id="svg-btn">SVG</button>
        <button class="download-btn" id="png-btn">PNG</button>
        <button class="download-btn" id="csv-btn">CSV</button>
        <button class="download-btn" id="excel-btn">Excel</button>
        <div class="export-error" style="display: none;"></div>
      </div>
    `
    vi.spyOn(document, `createElement`).mockImplementation((tag: string) =>
      tag === `a`
        ? ({ href: ``, download: ``, click: mock_click } as unknown as HTMLAnchorElement)
        : document.createElement(tag)
    )
  })

  afterEach(() => {
    vi.clearAllMocks()
    document.body.innerHTML = ``
  })

  it(`has all download buttons with correct structure and formatting`, () => {
    const buttons = document.querySelectorAll(`.download-btn`)
    expect(buttons).toHaveLength(4)

    const expected_formats = [`SVG`, `PNG`, `CSV`, `Excel`]
    buttons.forEach((button, idx) => {
      const format = expected_formats[idx]
      expect(button.textContent?.trim()).toBe(format)
      // Only check uppercase for SVG, PNG, CSV (not Excel)
      if (format !== `Excel`) {
        expect(button.textContent).toBe(button.textContent?.toUpperCase())
      }
    })
  })

  // Test all export formats with parameterized testing
  it.each([[`SVG`, `svg`], [`PNG`, `png`], [`CSV`, `csv`], [`Excel`, `excel`]] as const)(
    `triggers %s export correctly when button clicked`,
    async (format: string, id: string) => {
      const generate_spy = vi
        .fn()
        .mockResolvedValue({ filename: `test.${id}`, url: `mock-url` })
      const handle_export_fn = handle_export(generate_spy, format, { ...default_state })

      const button = document.querySelector(`#${id}-btn`) as HTMLButtonElement
      button.addEventListener(`click`, handle_export_fn)
      await button.click()

      expect(generate_spy).toHaveBeenCalledTimes(1)
      expect(generate_spy).toHaveBeenCalledWith({
        show_non_compliant: false,
        discovery_set: `unique_prototypes`,
      })

      await Promise.resolve() // Allow async handler to complete
      const error_el = document.querySelector(`.export-error`) as HTMLElement
      expect(error_el.style.display).toBe(`none`)
    },
  )

  // Test error handling scenarios with parameterized testing
  it.each([
    [
      `error`,
      new Error(`Failed to generate export`),
      `Error exporting Test: Failed to generate export`,
    ],
    [`null return`, null, `Failed to generate Test. The export function returned null.`],
  ])(`handles export %s correctly`, async (scenario, mock_result, expected_message) => {
    const mock_fn = scenario === `error`
      ? vi.fn().mockRejectedValue(mock_result)
      : vi.fn().mockResolvedValue(mock_result)

    const error_el = document.querySelector(`.export-error`) as HTMLElement
    const state = { ...default_state }

    const update_error_display = (err: string | null) => {
      state.export_error = err
      if (err) {
        error_el.textContent = err
        error_el.style.display = `block`
      } else {
        error_el.style.display = `none`
      }
    }

    const handle_fn = handle_export(mock_fn, `Test`, state)
    await handle_fn()

    if (state.export_error) {
      update_error_display(state.export_error)
    }

    expect(error_el.style.display).toBe(`block`)
    expect(error_el.textContent).toBe(expected_message)
  })

  it(`has error element hidden initially`, () => {
    const error_el = document.querySelector(`.export-error`) as HTMLElement
    expect(error_el.style.display).toBe(`none`)
  })
})
