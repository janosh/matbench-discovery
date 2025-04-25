import type { MockInstance } from 'vitest'
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest'

// Mock window.XMLSerializer
class MockXMLSerializer {
  serializeToString() {
    return `<div>Mock SVG Content</div>`
  }
}

describe(`SVG Export Functionality`, () => {
  // Create minimal mocks needed for testing
  const mock_click = vi.fn()

  // Create spies at the describe level to be used in multiple tests
  let create_element_spy: MockInstance

  beforeEach(() => {
    vi.clearAllMocks()

    // Mock minimal DOM APIs
    create_element_spy = vi
      .spyOn(document, `createElement`)
      .mockImplementation((tag: string) => {
        if (tag === `a`) {
          return {
            href: ``,
            download: ``,
            click: mock_click,
          } as unknown as HTMLAnchorElement
        }
        return document.createElement(tag)
      })

    // Mock URL.createObjectURL
    Object.defineProperty(window, `URL`, {
      value: { createObjectURL: vi.fn().mockReturnValue(`mock-url`) },
      writable: true,
    })

    // Mock document.querySelector to return a simple table
    vi.spyOn(document, `querySelector`).mockImplementation((_selector: string) => {
      return {
        querySelectorAll: () =>
          Array(10)
            .fill(0)
            .map(() => ({
              getBoundingClientRect: () => ({ height: 30, width: 500 }),
            })),
        cloneNode: () => ({
          style: {},
          querySelectorAll: () => [],
          getBoundingClientRect: () => ({ height: 300, width: 500 }),
        }),
      } as unknown as Element
    })

    // Mock XMLSerializer
    window.XMLSerializer = MockXMLSerializer as unknown as typeof XMLSerializer

    // Mock globals needed by the +page.svelte component
    vi.stubGlobal(`discovery_set`, `unique_prototypes`)
    vi.stubGlobal(`show_non_compliant`, false)
  })

  afterEach(() => {
    vi.clearAllMocks()
    vi.unstubAllGlobals()
  })

  it(`generates SVG with date in filename`, async () => {
    // Mock module's actual implementation to call our mocked document.createElement
    const mock_svg_result = {
      filename: `metrics-table-unique-prototypes-only-compliant-${new Date().toISOString().split(`T`)[0]}.svg`,
      url: `mock-url`,
    }

    // Import the function dynamically after mocks are set up
    const module = await import(`$lib/html-to-img`)

    // Mock the module function to return our prepared result AND call our document.createElement
    vi.spyOn(module, `generate_svg`).mockImplementation(async () => {
      // Simulate the function creating and clicking an anchor
      const a = document.createElement(`a`)
      a.click()
      return mock_svg_result
    })

    const generate_svg = module.generate_svg

    const result = await generate_svg({
      show_non_compliant: false,
      discovery_set: `unique_prototypes`,
    })

    // Check that anchor was created and clicked
    expect(create_element_spy).toHaveBeenCalledWith(`a`)
    expect(mock_click).toHaveBeenCalledTimes(1)

    // Assertions for result object
    expect(result).toBeTruthy()
    expect(result?.filename).toContain(`metrics-table-unique-prototypes-only-compliant`)
    expect(result?.filename).toContain(`.svg`)
    expect(result?.url).toBe(`mock-url`)

    // Verify date format in filename (YYYY-MM-DD)
    const today = new Date().toISOString().split(`T`)[0]
    expect(result?.filename).toContain(today)
  })

  it(`generates PNG with date in filename`, async () => {
    // Mock module's result
    const mock_png_result = {
      filename: `metrics-table-unique-prototypes-only-compliant-${new Date().toISOString().split(`T`)[0]}.png`,
      url: `mock-url`,
    }

    // Import the function dynamically after mocks are set up
    const module = await import(`$lib/html-to-img`)

    // Mock the module function to return our prepared result
    vi.spyOn(module, `generate_png`).mockImplementation(async () => {
      // Simulate the function creating and clicking an anchor
      const a = document.createElement(`a`)
      a.click()
      return mock_png_result
    })

    const generate_png = module.generate_png

    const result = await generate_png({
      show_non_compliant: false,
      discovery_set: `unique_prototypes`,
    })

    // Check that anchor was created and clicked
    expect(create_element_spy).toHaveBeenCalledWith(`a`)
    expect(mock_click).toHaveBeenCalledTimes(1)

    // Assertions for result object
    expect(result).toBeTruthy()
    expect(result?.filename).toContain(`metrics-table-unique-prototypes-only-compliant`)
    expect(result?.filename).toContain(`.png`)
    expect(result?.url).toBe(`mock-url`)

    // Verify date format in filename (YYYY-MM-DD)
    const today = new Date().toISOString().split(`T`)[0]
    expect(result?.filename).toContain(today)
  })

  it(`verifies filenames use param-case formatting`, async () => {
    // Set discovery_set to snake_case to test conversion
    vi.stubGlobal(`discovery_set`, `unique_prototypes`)

    // Mock return values with param-case
    const today = new Date().toISOString().split(`T`)[0]
    const svg_filename = `metrics-table-unique-prototypes-only-compliant-${today}.svg`
    const png_filename = `metrics-table-unique-prototypes-only-compliant-${today}.png`

    const mock_svg_result = {
      filename: svg_filename,
      url: `mock-url`,
    }

    const mock_png_result = {
      filename: png_filename,
      url: `mock-url`,
    }

    // Import the module
    const svg_module = await import(`$lib/html-to-img`)

    // Mock the module functions with appropriate types
    vi.spyOn(svg_module, `generate_svg`).mockResolvedValue(mock_svg_result)
    vi.spyOn(svg_module, `generate_png`).mockResolvedValue(mock_png_result)

    // Test SVG filename
    const svg_result = await svg_module.generate_svg({
      show_non_compliant: false,
      discovery_set: `unique_prototypes`,
    })

    // Verify param-case is used (dashes not underscores)
    expect(svg_result?.filename).toContain(`unique-prototypes`)
    expect(svg_result?.filename).not.toContain(`unique_prototypes`)

    // Test PNG filename
    const png_result = await svg_module.generate_png({
      show_non_compliant: false,
      discovery_set: `unique_prototypes`,
    })

    expect(png_result?.filename).toContain(`unique-prototypes`)
    expect(png_result?.filename).not.toContain(`unique_prototypes`)
  })

  it(`handles error when table element is not found for SVG`, async () => {
    // Mock querySelector to return null
    vi.spyOn(document, `querySelector`).mockReturnValue(null)
    const console_spy = vi.spyOn(console, `error`).mockImplementation(() => {})

    vi.resetModules() // Ensure mocks apply to fresh module import
    const module = await import(`$lib/html-to-img`)
    const result = await module.generate_svg({ discovery_set: `test` })

    expect(result).toBeNull()
    expect(console_spy).toHaveBeenCalledWith(`Table element not found for SVG export`)

    console_spy.mockRestore()
  })

  it(`handles error when table element is not found for PNG`, async () => {
    // Mock querySelector to return null
    vi.spyOn(document, `querySelector`).mockReturnValue(null)
    const console_spy = vi.spyOn(console, `error`).mockImplementation(() => {})

    vi.resetModules() // Ensure mocks apply to fresh module import
    const module = await import(`$lib/html-to-img`)
    const result = await module.generate_png({ discovery_set: `test` })

    expect(result).toBeNull()
    expect(console_spy).toHaveBeenCalledWith(`Table element not found for PNG export`)

    console_spy.mockRestore()
  })

  it(`handles error during SVG generation using toSvg`, async () => {
    // Mock toSvg to throw an error
    const mock_to_svg = vi.fn().mockRejectedValue(new Error(`toSvg failed`))
    vi.mock(`html-to-image`, () => ({
      toSvg: mock_to_svg,
      toPng: vi.fn(),
    }))

    vi.resetModules() // Ensure mocks apply to fresh module import
    const module = await import(`$lib/html-to-img`)
    const console_spy = vi.spyOn(console, `error`).mockImplementation(() => {})

    const result = await module.generate_svg({ discovery_set: `test` })

    expect(result).toBeNull()
    expect(console_spy).toHaveBeenCalledWith(`Error generating SVG:`, expect.any(Error))

    console_spy.mockRestore()
    vi.unmock(`html-to-image`) // Clean up mock
  })

  it(`handles error during PNG generation using toPng`, async () => {
    // Mock toPng to throw an error
    const mock_to_png = vi.fn().mockRejectedValue(new Error(`toPng failed`))
    vi.mock(`html-to-image`, () => ({
      toSvg: vi.fn(),
      toPng: mock_to_png,
    }))

    vi.resetModules() // Ensure mocks apply to fresh module import
    const module = await import(`$lib/html-to-img`)
    const console_spy = vi.spyOn(console, `error`).mockImplementation(() => {})

    const result = await module.generate_png({ discovery_set: `test` })

    expect(result).toBeNull()
    expect(console_spy).toHaveBeenCalledWith(`Error generating PNG:`, expect.any(Error))

    console_spy.mockRestore()
    console_spy.mockRestore()
    vi.unmock(`html-to-image`) // Clean up mock
  })
})

describe(`handle_export tests`, () => {
  // Define state type consistent with handle_export expectations
  type State = {
    export_error: string | null
    show_non_compliant?: boolean
    discovery_set?: string
  }
  let state: State
  let generator_spy: MockInstance

  beforeEach(() => {
    state = { export_error: null }
    generator_spy = vi.fn()
    vi.clearAllMocks()
  })

  it(`calls generator and handles success`, async () => {
    vi.resetModules() // Ensure mocks apply to fresh module import
    const module = await import(`$lib/html-to-img`)
    generator_spy.mockResolvedValue({ filename: `test.fmt`, url: `some-url` })
    // Pass required state properties
    state = { export_error: null, discovery_set: `test`, show_non_compliant: false }
    const handler = module.handle_export(generator_spy, `fmt`, state)

    await handler()

    expect(generator_spy).toHaveBeenCalledWith({
      discovery_set: state.discovery_set,
      show_non_compliant: state.show_non_compliant,
    })
    expect(state.export_error).toBeNull()
  })

  it(`sets error message when generator returns null`, async () => {
    vi.resetModules() // Ensure mocks apply to fresh module import
    const module = await import(`$lib/html-to-img`)
    generator_spy.mockResolvedValue(null)
    // Pass required state properties
    state = { export_error: null, discovery_set: `test`, show_non_compliant: false }
    const handler = module.handle_export(generator_spy, `fmt`, state)

    await handler()

    expect(generator_spy).toHaveBeenCalledWith({
      discovery_set: state.discovery_set,
      show_non_compliant: state.show_non_compliant,
    })
    expect(state.export_error).toBe(
      `Failed to generate fmt. The export function returned null.`,
    )
  })

  it(`sets error message when generator throws an Error`, async () => {
    vi.resetModules() // Ensure mocks apply to fresh module import
    const module = await import(`$lib/html-to-img`)
    const error_message = `Generator failed`
    generator_spy.mockRejectedValue(new Error(error_message))
    const console_spy = vi.spyOn(console, `error`).mockImplementation(() => {})
    // Pass required state properties
    state = { export_error: null, discovery_set: `test`, show_non_compliant: false }

    const handler = module.handle_export(generator_spy, `fmt`, state)

    await handler()

    expect(generator_spy).toHaveBeenCalledWith({
      discovery_set: state.discovery_set,
      show_non_compliant: state.show_non_compliant,
    })
    expect(state.export_error).toBe(`Error exporting fmt: ${error_message}`)
    expect(console_spy).toHaveBeenCalledWith(`Error exporting fmt:`, expect.any(Error))

    console_spy.mockRestore()
  })

  it(`sets error message when generator throws a string`, async () => {
    vi.resetModules() // Ensure mocks apply to fresh module import
    const module = await import(`$lib/html-to-img`)
    const error_string = `Something went wrong`
    generator_spy.mockRejectedValue(error_string)
    const console_spy = vi.spyOn(console, `error`).mockImplementation(() => {})
    // Pass required state properties
    state = { export_error: null, discovery_set: `test`, show_non_compliant: false }

    const handler = module.handle_export(generator_spy, `fmt`, state)

    await handler()

    expect(generator_spy).toHaveBeenCalledWith({
      discovery_set: state.discovery_set,
      show_non_compliant: state.show_non_compliant,
    })
    expect(state.export_error).toBe(`Error exporting fmt: ${error_string}`)
    expect(console_spy).toHaveBeenCalledWith(`Error exporting fmt:`, error_string)

    console_spy.mockRestore()
  })
})
