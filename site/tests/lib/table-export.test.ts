import {
  generate_csv,
  generate_excel,
  generate_png,
  generate_svg,
  handle_export,
} from '$lib/table-export'
import { toPng, toSvg } from 'html-to-image'
import { beforeEach, describe, expect, it, vi } from 'vitest'

vi.mock(`html-to-image`, () => ({
  toSvg: vi.fn(),
  toPng: vi.fn(),
}))

const today = () => new Date().toISOString().split(`T`)[0]

const mount_table = (): HTMLTableElement => {
  document.body.innerHTML = `
    <table class="heatmap">
      <thead>
        <tr>
          <th>Model</th>
          <th>F1</th>
          <th>Org</th>
          <th>DAF</th>
          <th>Links</th>
          <th>CPS</th>
          <th>R<sup>2</sup></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td data-sort-value="Model A"><a href="/models/a">Model A</a></td>
          <td data-sort-value="0.851">0.851</td>
          <td><svg aria-label="Org"></svg></td>
          <td data-sort-value="2.346">2.346</td>
          <td><a href="https://example.com">Link</a></td>
          <td data-sort-value="1.235">1.235</td>
          <td>0.75</td>
        </tr>
        <tr>
          <td data-sort-value='Model "Special"'>Model "Special"</td>
          <td data-sort-value="0.123">0.123</td>
          <td><svg aria-label="Org"></svg></td>
          <td data-sort-value="1.234">1.234</td>
          <td><a href="https://example.org">Link</a></td>
          <td data-sort-value="0.987">0.987</td>
          <td>0.25</td>
        </tr>
      </tbody>
    </table>
  `

  const table = document.querySelector<HTMLTableElement>(`.heatmap`)
  if (!table) throw new Error(`Failed to mount test table`)
  return table
}

describe(`Table Export Functionality`, () => {
  let click_spy: ReturnType<typeof vi.spyOn>
  let create_object_url_spy: ReturnType<typeof vi.spyOn>

  beforeEach(() => {
    vi.restoreAllMocks()
    mount_table()

    vi.mocked(toSvg).mockResolvedValue(`data:image/svg+xml;base64,test`)
    vi.mocked(toPng).mockResolvedValue(`data:image/png;base64,test`)
    click_spy = vi
      .spyOn(HTMLAnchorElement.prototype, `click`)
      .mockImplementation(() => {})
    globalThis.URL.createObjectURL = (() => `mock-url`) as typeof URL.createObjectURL
    globalThis.URL.revokeObjectURL = (() => {}) as typeof URL.revokeObjectURL
    create_object_url_spy = vi
      .spyOn(globalThis.URL, `createObjectURL`)
      .mockReturnValue(`mock-url`)
  })

  it.each([
    [`SVG`, generate_svg, toSvg, `.svg`],
    [`PNG`, generate_png, toPng, `.png`],
  ] as const)(
    `generates %s export with a cleaned table clone`,
    async (_format, generator, encoder, extension) => {
      let captured_container: HTMLElement | undefined
      vi.mocked(encoder).mockImplementation((container) => {
        captured_container = container
        return Promise.resolve(`data:image/${extension.slice(1)};base64,test`)
      })

      const result = await generator({
        show_non_compliant: false,
        discovery_set: `unique_prototypes`,
      })

      if (!result) throw new Error(`${_format} export returned null`)
      if (!captured_container) throw new Error(`${_format} export did not call encoder`)

      expect(result.filename).toBe(
        `matbench-unique-prototypes-compliant-2models-${today()}${extension}`,
      )
      expect(result.url).toBe(extension === `.svg` ? `mock-url` : `data:image/png;base64,test`)
      expect(click_spy).toHaveBeenCalled()
      expect(captured_container.querySelectorAll(`svg, a[href]`)).toHaveLength(0)
      expect(captured_container.querySelectorAll(`sub, sup`)).toHaveLength(1)
    },
  )

  it.each([
    [`SVG`, generate_svg],
    [`PNG`, generate_png],
    [`CSV`, generate_csv],
    [`Excel`, generate_excel],
  ] as const)(`returns null when table is missing for %s export`, async (_format, generator) => {
    document.body.innerHTML = ``
    const console_spy = vi.spyOn(console, `error`).mockImplementation(() => {})

    await expect(Promise.resolve(generator({ discovery_set: `test` }))).resolves.toBeNull()
    expect(console_spy).toHaveBeenCalled()
  })

  it(`generates CSV with formatted data and excludes icon columns`, async () => {
    const result = generate_csv({
      show_non_compliant: false,
      discovery_set: `unique_prototypes`,
    })

    if (!result) throw new Error(`CSV export returned null`)

    expect(result.filename).toBe(`matbench-unique-prototypes-compliant-2models-${today()}.csv`)
    expect(result.url).toBe(`mock-url`)
    expect(click_spy).toHaveBeenCalled()

    const blob = create_object_url_spy.mock.calls[0]?.[0]
    expect(blob).toBeInstanceOf(Blob)
    const csv_content = await (blob as Blob).text()
    expect(csv_content).toContain(`Model,F1,DAF,CPS,R2`)
    expect(csv_content).not.toContain(`Org`)
    expect(csv_content).not.toContain(`Links`)
    expect(csv_content).toContain(`Model A`)
    expect(csv_content).toContain(`"Model ""Special"""`)
  })

  it(`generates Excel with the expected MIME type`, async () => {
    const result = await generate_excel({
      show_non_compliant: true,
      discovery_set: `test_set`,
    })

    if (!result) throw new Error(`Excel export returned null`)

    expect(result.filename).toBe(`matbench-test-set-all-2models-${today()}.xlsx`)
    expect(result.url).toBe(`mock-url`)
    expect(click_spy).toHaveBeenCalled()

    const blob = create_object_url_spy.mock.calls[0]?.[0]
    expect(blob).toBeInstanceOf(Blob)
    expect((blob as Blob).type).toBe(
      `application/vnd.openxmlformats-officedocument.spreadsheetml.sheet`,
    )
  })

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
        const handler = handle_export(generator_spy, `fmt`, state)
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
