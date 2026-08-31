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
  // includes the structural rank (#) column MetricsTable renders via
  // show_row_numbers, which exports must skip like the icon columns
  document.body.innerHTML = `
    <table class="heatmap">
      <thead>
        <tr>
          <th class="row-num-col">#</th>
          <th>Model</th>
          <th>F1</th>
          <th>Org</th>
          <th>DAF</th>
          <th>Links</th>
          <th>CPS</th>
          <th>R<sup>2</sup></th>
          <th>Params</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td class="row-num-col">1</td>
          <td data-sort-value="Model A"><a href="/models/a">Model A</a></td>
          <td data-sort-value="0.851">0.851</td>
          <td><svg aria-label="Org"></svg></td>
          <td data-sort-value="2.346">2.346</td>
          <td><a href="https://example.com">Link</a></td>
          <td data-sort-value="1.235">1.235</td>
          <td>0.75</td>
          <td><span title="4,700,000 params" data-sort-value="4700000">4.7M</span></td>
        </tr>
        <tr>
          <td class="row-num-col">2</td>
          <td data-sort-value='Model "Special"'>Model "Special"</td>
          <td data-sort-value="0.123">0.123</td>
          <td><svg aria-label="Org"></svg></td>
          <td data-sort-value="1.234">1.234</td>
          <td><a href="https://example.org">Link</a></td>
          <td data-sort-value="0.987">0.987</td>
          <td>0.25</td>
          <td><span data-sort-value="900">900</span></td>
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
    mount_table()

    vi.mocked(toSvg).mockResolvedValue(`data:image/svg+xml;base64,test`)
    vi.mocked(toPng).mockResolvedValue(`data:image/png;base64,test`)
    click_spy = vi
      .spyOn(HTMLAnchorElement.prototype, `click`)
      .mockImplementation(() => {})
    globalThis.URL.createObjectURL = () => `mock-url`
    globalThis.URL.revokeObjectURL = () => {}
    create_object_url_spy = vi
      .spyOn(globalThis.URL, `createObjectURL`)
      .mockReturnValue(`mock-url`)
  })

  // Blob handed to URL.createObjectURL during an export (i.e. the downloaded file contents)
  const exported_blob = (): Blob => {
    const blob = create_object_url_spy.mock.calls[0]?.[0]
    expect(blob).toBeInstanceOf(Blob)
    return blob as Blob
  }

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

      const result = await generator({ discovery_set: `unique_prototypes` })

      if (!result) throw new Error(`${_format} export returned null`)
      if (!captured_container) throw new Error(`${_format} export did not call encoder`)

      expect(result.filename).toBe(
        `matbench-unique-prototypes-2models-${today()}${extension}`,
      )
      expect(result.url).toBe(
        extension === `.svg` ? `mock-url` : `data:image/png;base64,test`,
      )
      expect(click_spy).toHaveBeenCalled()
      expect(captured_container.querySelectorAll(`svg, a[href]`)).toHaveLength(0)
      expect(captured_container.querySelectorAll(`sub, sup`)).toHaveLength(1)
      // the structural rank (#) column is stripped from image exports
      expect(captured_container.querySelector(`.row-num-col`)).toBeNull()
    },
  )

  it.each([
    [`SVG`, generate_svg],
    [`PNG`, generate_png],
    [`CSV`, generate_csv],
    [`Excel`, generate_excel],
  ] as const)(
    `returns null when table is missing for %s export`,
    async (_format, generator) => {
      document.body.innerHTML = ``
      const console_spy = vi.spyOn(console, `error`).mockImplementation(() => {})

      await expect(
        Promise.resolve(generator({ discovery_set: `test` })),
      ).resolves.toBeNull()
      expect(console_spy).toHaveBeenCalled()
    },
  )

  it(`generates CSV with formatted data and excludes icon columns`, async () => {
    const result = generate_csv({ discovery_set: `unique_prototypes` })

    if (!result) throw new Error(`CSV export returned null`)

    expect(result.filename).toBe(`matbench-unique-prototypes-2models-${today()}.csv`)
    expect(result.url).toBe(`mock-url`)
    expect(click_spy).toHaveBeenCalled()

    const csv_content = await exported_blob().text()
    expect(csv_content).toContain(`Model,F1,DAF,CPS,R2,Params`)
    expect(csv_content).not.toContain(`Org`)
    expect(csv_content).not.toContain(`Links`)
    expect(csv_content).not.toContain(`#`) // rank column excluded
    expect(csv_content).toContain(`Model A`)
    expect(csv_content).toContain(`"Model ""Special"""`)
    // HTML-string cells carry their raw value on an inner span, not the <td>: export the
    // unformatted number (rendered per the Params '~s' spec) rather than the display text
    expect(csv_content).toContain(`4.70e+6`)
    expect(csv_content).not.toContain(`4.7M`)
    expect(csv_content).toContain(`,900`)
  })

  it(`cleans text cells: ignores data-sort-value="null", collapses whitespace, decodes entities`, async () => {
    document.body.innerHTML = `
      <table class="heatmap">
        <thead><tr><th>Model</th><th>Notes</th></tr></thead>
        <tbody>
          <tr>
            <td data-sort-value="null">spaced    out    text</td>
            <td>&amp;lt;tag&amp;gt; &amp;amp; more</td>
          </tr>
        </tbody>
      </table>
    `
    const result = generate_csv({ discovery_set: `test` })
    if (!result) throw new Error(`CSV export returned null`)

    const csv_content = await exported_blob().text()
    expect(csv_content).toContain(`Model,Notes`)
    // data-sort-value="null" is ignored, falling through to text with collapsed whitespace
    expect(csv_content).toContain(`spaced out text`)
    // double-encoded entities (&amp;lt; -> &lt; -> <) are decoded back to literals
    expect(csv_content).toContain(`<tag> & more`)
  })

  it(`generates Excel with the expected MIME type`, async () => {
    const result = await generate_excel({ discovery_set: `test_set` })

    if (!result) throw new Error(`Excel export returned null`)

    expect(result.filename).toBe(`matbench-test-set-2models-${today()}.xlsx`)
    expect(result.url).toBe(`mock-url`)
    expect(click_spy).toHaveBeenCalled()

    expect(exported_blob().type).toBe(
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
        const opts = { discovery_set: `test` }

        if (generator_result instanceof Error || typeof generator_result === `string`) {
          generator_spy.mockRejectedValue(generator_result)
        } else {
          generator_spy.mockResolvedValue(generator_result)
        }

        vi.spyOn(console, `error`).mockImplementation(() => {})
        expect(await handle_export(generator_spy, `fmt`, opts)).toBe(expected_error)
        expect(generator_spy).toHaveBeenCalledWith(opts)
      },
    )
  })
})
