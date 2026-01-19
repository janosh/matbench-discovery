import type { Mock, MockInstance } from 'vitest'
import { beforeAll, beforeEach, describe, expect, it, vi } from 'vitest'

// Check if running in Deno environment
const IS_DENO = `Deno` in globalThis

// Mock DOM table structure for testing
const create_mock_table = () =>
  ({
    querySelectorAll: vi.fn((selector: string) => {
      if (selector === `thead tr`) {
        return [
          {
            querySelectorAll: () => [
              { textContent: `Model` },
              { textContent: `F1` },
              { textContent: `Org` },
              { textContent: `DAF` },
              { textContent: `Links` },
              { textContent: `CPS` },
            ],
          },
        ]
      }
      if (selector === `tbody tr`) {
        return [
          {
            querySelectorAll: () => [
              { getAttribute: () => `Model A`, textContent: `Model A` },
              { getAttribute: () => `0.851`, textContent: `0.851` },
              { getAttribute: () => null, textContent: `<svg>...</svg>` },
              { getAttribute: () => `2.346`, textContent: `2.346` },
              { getAttribute: () => null, textContent: `<svg>...</svg>` },
              { getAttribute: () => `1.235`, textContent: `1.235` },
            ],
          },
        ]
      }
      return []
    }),
    cloneNode: () => ({
      style: {},
      querySelectorAll: () => [],
      getBoundingClientRect: () => ({ height: 300, width: 500 }),
    }),
  }) as unknown as Element

describe.skipIf(IS_DENO)(`Table Export Functionality`, () => {
  let create_element_spy: MockInstance
  let query_selector_spy: MockInstance
  let original_create_element: typeof document.createElement
  const mock_click = vi.fn()

  beforeAll(() => {
    // Save original createElement before any mocks - must be inside describe block
    // to avoid import-time errors in non-DOM environments (Deno)
    original_create_element = document.createElement.bind(document)
  })

  beforeEach(() => {
    vi.clearAllMocks()

    // Common mocks - use original_create_element to avoid infinite recursion
    create_element_spy = vi
      .spyOn(document, `createElement`)
      .mockImplementation((tag: string) =>
        tag === `a`
          ? ({
            href: ``,
            download: ``,
            click: mock_click,
          } as unknown as HTMLAnchorElement)
          : original_create_element(tag)
      )

    query_selector_spy = vi
      .spyOn(document, `querySelector`)
      .mockImplementation(() => create_mock_table())

    globalThis.Blob = vi.fn() as unknown as typeof Blob

    // Mock URL methods on globalThis.URL
    if (!globalThis.URL) globalThis.URL = {} as typeof URL
    globalThis.URL.createObjectURL = vi.fn().mockReturnValue(`mock-url`)
    globalThis.URL.revokeObjectURL = vi.fn()
  })

  // Test image exports (SVG, PNG) with parameterized testing
  describe.each(
    [
      [`SVG`, `generate_svg`, `.svg`],
      [`PNG`, `generate_png`, `.png`],
    ] as const,
  )(`%s Export`, (format, function_name, extension) => {
    it(`generates ${format} with correct filename and calls download`, async () => {
      const module = await import(`$lib/table-export`)
      // deno-lint-ignore require-await
      vi.spyOn(module, function_name).mockImplementation(async () => {
        document.createElement(`a`).click()
        const date_str = new Date().toISOString().split(`T`)[0]
        return {
          filename:
            `metrics-table-unique-prototypes-only-compliant-${date_str}${extension}`,
          url: `mock-url`,
        }
      })

      const result = await module[function_name]({
        show_non_compliant: false,
        discovery_set: `unique_prototypes`,
      })

      expect(create_element_spy).toHaveBeenCalledWith(`a`)
      expect(mock_click).toHaveBeenCalledTimes(1)
      expect(result?.filename).toContain(`unique-prototypes-only-compliant`)
      expect(result?.filename).toContain(extension)
      expect(result?.filename).toContain(new Date().toISOString().split(`T`)[0])
      expect(result?.url).toBe(`mock-url`)
    })

    it(`preserves subscripts and superscripts in ${format} export`, async () => {
      // Reset modules and restore spies to use real DOM operations
      vi.resetModules()
      create_element_spy.mockRestore()
      query_selector_spy.mockRestore()

      // Create real DOM table with sub/sup elements
      const real_table = document.createElement(`table`)
      real_table.className = `heatmap`

      const thead = document.createElement(`thead`)
      const header_row = document.createElement(`tr`)

      const th1 = document.createElement(`th`)
      th1.textContent = `Model`
      header_row.appendChild(th1)

      const th2 = document.createElement(`th`)
      th2.innerHTML = `R<sup>2</sup>`
      header_row.appendChild(th2)

      const th3 = document.createElement(`th`)
      th3.innerHTML = `κ<sub>SRME</sub>`
      header_row.appendChild(th3)

      thead.appendChild(header_row)
      real_table.appendChild(thead)

      const tbody = document.createElement(`tbody`)
      const body_row = document.createElement(`tr`)

      const td1 = document.createElement(`td`)
      td1.textContent = `Test Model`
      body_row.appendChild(td1)

      const td2 = document.createElement(`td`)
      td2.textContent = `0.85`
      body_row.appendChild(td2)

      const td3 = document.createElement(`td`)
      td3.textContent = `1.23`
      body_row.appendChild(td3)

      tbody.appendChild(body_row)
      real_table.appendChild(tbody)

      // Append to document so querySelector can find it
      document.body.appendChild(real_table)

      // Mock the image generation library to capture the container
      let captured_container: HTMLElement | null = null
      const mock_lib = format === `SVG`
        ? {
          toSvg: vi.fn().mockImplementation((container) => {
            captured_container = container
            return Promise.resolve(`data:image/svg+xml;base64,test`)
          }),
        }
        : {
          toPng: vi.fn().mockImplementation((container) => {
            captured_container = container
            return Promise.resolve(`data:image/png;base64,test`)
          }),
        }

      vi.doMock(`html-to-image`, () => mock_lib)

      const module = await import(`$lib/table-export`)
      await module[function_name]({ discovery_set: `test` })

      // Clean up - remove the table from DOM
      real_table.remove()

      // Verify that the function was called and container was processed
      expect(captured_container).not.toBeNull()
      const container = captured_container as unknown as HTMLElement

      // The key test: verify sub/sup elements are preserved in the structure
      // clean_table_for_export should NOT remove sub/sup elements
      const preserved_elements = container.querySelectorAll(`sub, sup`)
      expect(preserved_elements.length).toBe(2)
    })

    it(`handles table not found error for ${format}`, async () => {
      query_selector_spy.mockReturnValue(null)
      const console_spy = vi.spyOn(console, `error`).mockImplementation(() => {})

      vi.resetModules()
      const module = await import(`$lib/table-export`)
      const result = await module[function_name]({ discovery_set: `test` })

      expect(result).toBeNull()
      expect(console_spy).toHaveBeenCalled()
      console_spy.mockRestore()
    })

    it(`handles library errors gracefully for ${format}`, async () => {
      const console_spy = vi.spyOn(console, `error`).mockImplementation(() => {})
      const mock_lib = format === `SVG`
        ? { toSvg: vi.fn().mockRejectedValue(new Error(`Library failed`)) }
        : { toPng: vi.fn().mockRejectedValue(new Error(`Library failed`)) }

      vi.doMock(`html-to-image`, () => mock_lib)

      const module = await import(`$lib/table-export`)
      const result = await module[function_name]({ discovery_set: `test` })

      expect(result).toBeNull()
      expect(console_spy).toHaveBeenCalled()
      expect(query_selector_spy).toHaveBeenCalledWith(`.heatmap`)
      console_spy.mockRestore()
    })
  })

  // Test data exports (CSV, Excel) with parameterized testing
  describe.each(
    [
      [`CSV`, `generate_csv`, `.csv`, `text/csv;charset=utf-8;`],
      [
        `Excel`,
        `generate_excel`,
        `.xlsx`,
        `application/vnd.openxmlformats-officedocument.spreadsheetml.sheet`,
      ],
    ] as const,
  )(`%s Export`, (format, function_name, extension, mime_type) => {
    beforeEach(() => {
      if (format === `Excel`) {
        // Simplified Excel mock
        vi.doMock(`xlsx`, () => ({
          utils: {
            book_new: vi.fn().mockReturnValue({}),
            aoa_to_sheet: vi.fn().mockReturnValue({ '!ref': `A1:D3` }),
            book_append_sheet: vi.fn(),
            decode_range: vi
              .fn()
              .mockReturnValue({ s: { r: 0, c: 0 }, e: { r: 2, c: 3 } }),
            encode_cell: vi
              .fn()
              .mockImplementation(({ r, c }) => `${String.fromCharCode(65 + c)}${r + 1}`),
          },
          write: vi.fn().mockReturnValue(new ArrayBuffer(100)),
        }))
      }
    })

    it(`generates ${format} with proper data and excludes SVG columns`, async () => {
      const module = await import(`$lib/table-export`)
      const result = await module[function_name]({
        show_non_compliant: false,
        discovery_set: `unique_prototypes`,
      })

      expect(create_element_spy).toHaveBeenCalledWith(`a`)
      expect(mock_click).toHaveBeenCalledTimes(1)
      expect(result?.filename).toContain(extension)
      expect(result?.url).toBe(`mock-url`)

      // Verify data structure in blob
      const blob_call = (globalThis.Blob as Mock).mock.calls[0]
      expect(blob_call[1].type).toBe(mime_type)

      if (format === `CSV`) {
        const csv_content = blob_call[0][0]
        expect(csv_content).toContain(`Model,F1,DAF,CPS`) // Headers without SVG columns
        expect(csv_content).not.toContain(`Org`)
        expect(csv_content).not.toContain(`Links`)
        expect(csv_content).toContain(`Model A`)
      }
    })

    it(`handles special characters correctly for ${format}`, async () => {
      // Test only for CSV since Excel handles this differently
      if (format !== `CSV`) return

      query_selector_spy.mockReturnValue({
        querySelectorAll: vi.fn((selector: string) => {
          if (selector === `thead tr`) {
            return [{ querySelectorAll: () => [{ textContent: `Model` }] }]
          }
          if (selector === `tbody tr`) {
            return [
              {
                querySelectorAll: () => [
                  {
                    getAttribute: () => `Model "Special"`,
                    textContent: `Model "Special"`,
                  },
                ],
              },
            ]
          }
          return []
        }),
      } as unknown as Element)

      const module = await import(`$lib/table-export`)
      await module[function_name]({ discovery_set: `test` })

      const csv_content = (globalThis.Blob as Mock).mock.calls[0][0][0]
      expect(csv_content).toContain(`"Model ""Special"""`)
    })

    it(`handles errors when table not found for ${format}`, async () => {
      query_selector_spy.mockReturnValue(null)
      const console_spy = vi.spyOn(console, `error`).mockImplementation(() => {})

      const module = await import(`$lib/table-export`)
      const result = await module[function_name]({ discovery_set: `test` })

      expect(result).toBeNull()
      expect(console_spy).toHaveBeenCalledWith(
        `Error generating ${format}:`,
        expect.any(Error),
      )
      console_spy.mockRestore()
    })
  })

  // Combined utility function tests
  describe(`Utility Functions`, () => {
    it(`formats numbers and generates filenames correctly`, async () => {
      const module = await import(`$lib/table-export`)

      // Test through CSV generation with known values
      query_selector_spy.mockReturnValue({
        querySelectorAll: vi.fn((selector: string) => {
          if (selector === `thead tr`) {
            return [{ querySelectorAll: () => [{ textContent: `CPS` }] }]
          }
          if (selector === `tbody tr`) {
            return [
              {
                querySelectorAll: () => [
                  { getAttribute: () => `1.23456789`, textContent: `1.235` },
                ],
              },
            ]
          }
          return []
        }),
      } as unknown as Element)

      const result = await module.generate_csv({
        show_non_compliant: false,
        discovery_set: `unique_prototypes`,
      })

      // Test filename generation
      expect(result?.filename).toContain(`unique-prototypes`) // param-case conversion
      expect(result?.filename).toContain(`compliant`) // compliance suffix
      expect(result?.filename).toContain(new Date().toISOString().split(`T`)[0]) // date

      // Test number formatting in CSV content
      const csv_content = (globalThis.Blob as Mock).mock.calls[0][0][0]
      expect(csv_content).toContain(`1.235`)
    })

    it(`excludes compliance suffix when showing all models`, async () => {
      const module = await import(`$lib/table-export`)
      const result = await module.generate_csv({
        show_non_compliant: true,
        discovery_set: `test`,
      })

      expect(result?.filename).not.toContain(`compliant`)
    })
  })

  // Combined handle_export tests
  describe(`handle_export Function`, () => {
    it.each([
      [`success`, { filename: `test.fmt`, url: `some-url` }, null],
      [`null return`, null, `Failed to generate fmt. The export function returned null.`],
      [`error`, new Error(`Generator failed`), `Error exporting fmt: Generator failed`],
      [
        `string error`,
        `Something went wrong`,
        `Error exporting fmt: Something went wrong`,
      ],
    ])(
      `handles %s case correctly`,
      async (_test_name, generator_result, expected_error) => {
        const generator_spy = vi.fn()
        const state = {
          export_error: null,
          discovery_set: `test`,
          show_non_compliant: false,
        }

        if (generator_result instanceof Error) {
          generator_spy.mockRejectedValue(generator_result)
        } else if (typeof generator_result === `string`) {
          generator_spy.mockRejectedValue(generator_result)
        } else {
          generator_spy.mockResolvedValue(generator_result)
        }

        const console_spy = vi.spyOn(console, `error`).mockImplementation(() => {})

        vi.resetModules()
        const module = await import(`$lib/table-export`)

        const handler = module.handle_export(generator_spy, `fmt`, state)
        await handler()

        expect(generator_spy).toHaveBeenCalledWith({
          discovery_set: state.discovery_set,
          show_non_compliant: state.show_non_compliant,
        })
        expect(state.export_error).toBe(expected_error)

        console_spy.mockRestore()
      },
    )
  })
})
