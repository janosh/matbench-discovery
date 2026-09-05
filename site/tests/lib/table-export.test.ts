import { generate_csv } from '$lib/table-export'
import { beforeEach, describe, expect, it, vi } from 'vitest'

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

  it(`returns null when the table is missing`, () => {
    document.body.innerHTML = ``
    const console_spy = vi.spyOn(console, `error`).mockImplementation(() => {})
    expect(generate_csv({ discovery_set: `test` })).toBeNull()
    expect(console_spy).toHaveBeenCalled()
  })

  it(`generates CSV with formatted data and excludes icon columns`, async () => {
    const model_cell = document.querySelector<HTMLElement>(`tbody td:nth-child(2)`)
    if (!model_cell) throw new Error(`missing model cell`)
    model_cell.dataset.sortValue = `Model A, revised`
    vi.useFakeTimers()
    const revoke_spy = vi.spyOn(URL, `revokeObjectURL`)
    const result = generate_csv({ discovery_set: `unique_prototypes` })

    if (!result) throw new Error(`CSV export returned null`)

    expect(result.filename).toBe(`matbench-unique-prototypes-2models-${today()}.csv`)
    expect(result.url).toBe(`mock-url`)
    expect(click_spy).toHaveBeenCalled()
    expect(revoke_spy).not.toHaveBeenCalled()
    vi.advanceTimersByTime(100)
    expect(revoke_spy).toHaveBeenCalledWith(result.url)
    vi.useRealTimers()

    const csv_content = await exported_blob().text()
    expect(csv_content).toContain(`Model,F1,DAF,CPS,R2,Params`)
    expect(csv_content).not.toContain(`Org`)
    expect(csv_content).not.toContain(`Links`)
    expect(csv_content).not.toContain(`#`) // rank column excluded
    expect(csv_content).toContain(`"Model A, revised"`)
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
})
