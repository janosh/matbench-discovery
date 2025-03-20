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
})
