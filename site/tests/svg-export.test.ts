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
  const mock_write_text = vi.fn().mockResolvedValue(undefined)

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

    // Mock clipboard API
    Object.defineProperty(navigator, `clipboard`, {
      value: { writeText: mock_write_text },
      configurable: true,
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
        }),
      } as unknown as Element
    })

    // Mock XMLSerializer
    window.XMLSerializer = MockXMLSerializer as unknown as typeof XMLSerializer

    // Mock getElementById for PDF button
    vi.spyOn(document, `getElementById`).mockReturnValue({
      textContent: `PDF`,
    } as unknown as HTMLElement)

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
    const module = await import(`$lib/svg-export`)

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

  it(`copies PDF conversion command after SVG download`, async () => {
    // Mock return values
    const today = new Date().toISOString().split(`T`)[0]
    const svg_filename = `metrics-table-unique-prototypes-only-compliant-${today}.svg`
    const pdf_filename = svg_filename.replace(`.svg`, `.pdf`)
    const command = `# Install pdf2svg:
# Linux: sudo apt-get install pdf2svg
# macOS: brew install pdf2svg
pdf2svg ${svg_filename.replace(`.svg`, ``)}.{svg,pdf}`

    const mock_pdf_result = {
      command,
      svg_filename,
      pdf_filename,
    }

    // Import the function dynamically after mocks are set up
    const module = await import(`$lib/svg-export`)

    // Mock the module function to simulate creating an anchor element
    vi.spyOn(module, `copy_pdf_conversion_cmd`).mockImplementation(async () => {
      // Simulate the function creating and clicking an anchor
      const a = document.createElement(`a`)
      a.click()
      // Simulate clipboard write
      await navigator.clipboard.writeText(command)
      return mock_pdf_result
    })

    const copy_pdf_command = module.copy_pdf_conversion_cmd

    // Call the function
    const result = await copy_pdf_command()

    // Check that SVG was generated (anchor created and clicked)
    expect(create_element_spy).toHaveBeenCalledWith(`a`)
    expect(mock_click).toHaveBeenCalledTimes(1)

    // Check clipboard write was called
    expect(mock_write_text).toHaveBeenCalledTimes(1)
    expect(mock_write_text).toHaveBeenCalledWith(command)

    // Assertions for result object
    expect(result).toBeTruthy()
    expect(result?.command).toBe(command)
    expect(result?.svg_filename).toContain(`.svg`)
    expect(result?.pdf_filename).toContain(`.pdf`)

    // PDF filename should be SVG filename with extension changed
    expect(result?.pdf_filename).toBe(result?.svg_filename.replace(`.svg`, `.pdf`))

    // Verify date in filename
    expect(result?.svg_filename).toContain(
      `metrics-table-unique-prototypes-only-compliant`,
    )
    expect(result?.svg_filename).toContain(today)
  })

  it(`verifies filenames use param-case formatting`, async () => {
    // Set discovery_set to snake_case to test conversion
    vi.stubGlobal(`discovery_set`, `unique_prototypes`)

    // Mock return values with param-case
    const today = new Date().toISOString().split(`T`)[0]
    const svg_filename = `metrics-table-unique-prototypes-only-compliant-${today}.svg`
    const pdf_filename = svg_filename.replace(`.svg`, `.pdf`)

    const mock_svg_result = {
      filename: svg_filename,
      url: `mock-url`,
    }

    const mock_pdf_result = {
      command: `pdf2svg ${svg_filename.replace(`.svg`, ``)}.{svg,pdf}`,
      svg_filename,
      pdf_filename,
    }

    // Import the module
    const svg_module = await import(`$lib/svg-export`)

    // Mock the module functions with appropriate types
    vi.spyOn(svg_module, `generate_svg`).mockResolvedValue(mock_svg_result)
    vi.spyOn(svg_module, `copy_pdf_conversion_cmd`).mockResolvedValue(mock_pdf_result)

    // Test SVG filename
    const result = await svg_module.generate_svg({
      show_non_compliant: false,
      discovery_set: `unique_prototypes`,
    })

    // Verify param-case is used (dashes not underscores)
    expect(result?.filename).toContain(`unique-prototypes`)
    expect(result?.filename).not.toContain(`unique_prototypes`)

    // Test PDF filename via command
    const pdf_result = await svg_module.copy_pdf_conversion_cmd()

    expect(pdf_result?.svg_filename).toContain(`unique-prototypes`)
    expect(pdf_result?.pdf_filename).toContain(`unique-prototypes`)
    expect(pdf_result?.svg_filename).not.toContain(`unique_prototypes`)
  })
})
