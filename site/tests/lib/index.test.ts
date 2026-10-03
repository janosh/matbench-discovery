import data_files from '$pkg/data-files.yml'
import DATASETS from '$data/datasets.yml'
import { arr_to_str, format_date } from '#lib'
import { MODELS } from '#lib/models.svelte.js'
import { scatter_options_by_key } from '#lib/labels.js'
import { render_data_markdown } from '../../scripts/markdown-data.js'
import { sync_url_params } from '#lib/url-state.svelte.js'
import { goto } from '$app/navigation'
import { valid_query_param } from 'svelte-widgets/url-params'
import { describe, expect, it, vi } from 'vite-plus/test'

describe(`#lib data includes rendered YAML Markdown`, () => {
  it(`DATASETS entries expose computed slug and description_html`, () => {
    const entries = Object.entries(DATASETS)
    expect(entries.length).toBeGreaterThanOrEqual(20) // datasets.yml entry count
    for (const [key, dataset] of entries) {
      for (const field of [`slug`, `description_html`] as const)
        expect(dataset[field]?.length, `${key} missing ${field}`).toBeGreaterThan(0)
      expect(Object.keys(dataset.notes_html ?? {})).toEqual(
        Object.keys(dataset.notes ?? {}),
      )
    }
    expect(DATASETS.ELEMENTA.notes_html?.Access).toContain(
      `href="https://huggingface.co/datasets/kairosmaterial/ELEMENTA"`,
    )
  })

  it(`data_files entries expose computed html`, () => {
    const entries = Object.entries(data_files).filter(
      ([key, entry]) => !key.startsWith(`_`) && typeof entry === `object`,
    )
    expect(entries.length).toBeGreaterThanOrEqual(10) // data-files.yml entry count
    for (const [key, entry] of entries) {
      expect(
        (entry as { html?: string }).html?.length,
        `${key} missing html`,
      ).toBeGreaterThan(0)
    }
    const entry = data_files.wbm_computed_structure_entries
    expect(typeof entry === `object` && entry.html).toContain(
      `href="https://github.com/materialsproject/pymatgen/`,
    )
  })

  it(`model notes are rendered through eager YAML imports`, () => {
    const notes = MODELS.flatMap((model) => (model.notes ? [model.notes] : []))
    expect(notes.length).toBeGreaterThan(20)
    for (const note of notes) {
      for (const [key, value] of Object.entries(note)) {
        if (typeof value === `string` && value.trim())
          expect(
            note.html?.[key]?.length,
            `missing rendered note ${key}`,
          ).toBeGreaterThan(0)
      }
    }
  })

  it(`renders references while dropping raw HTML and retaining authored HTML overrides`, async () => {
    const files = {
      example: { description: `Before <b>bold</b> and [reference][target]` },
      _links: `[target]: https://example.org/reference`,
      _private: { description: `Leave **metadata** alone` },
    }
    await render_data_markdown(files, `/repo/matbench_discovery/data-files.yml`)
    expect(files.example).toEqual({
      description: `Before <b>bold</b> and [reference][target]`,
      html: `<p>Before bold and <a href="https://example.org/reference">reference</a></p>\n`,
    })
    expect(files._private).toEqual({ description: `Leave **metadata** alone` })
    const model = {
      notes: {
        description: `<div>Hidden</div>\n\n**Visible** {value}`,
        training: `Original`,
        constructor: `**Constructor**`,
        [`__proto__`]: `*Prototype*`,
        count: 2,
        html: { training: `<p>Authored override</p>` },
      },
    }
    await render_data_markdown(model, `/repo/models/example/model.yml`)
    expect(model.notes.html).toEqual({
      description: `<p><strong>Visible</strong> &#123;value&#125;</p>\n`,
      training: `<p>Authored override</p>`,
      constructor: `<p><strong>Constructor</strong></p>\n`,
      [`__proto__`]: `<p><em>Prototype</em></p>\n`,
    })
    expect(Object.getPrototypeOf(model.notes.html)).toBe(Object.prototype)
    const unrelated = { description: `**untouched**` }
    expect(await render_data_markdown(unrelated, `/repo/data/other.yml`)).toBe(unrelated)
    expect(unrelated).toEqual({ description: `**untouched**` })
  })

  it.each([
    [`---\n\n**Visible**`, `<hr>\n<p><strong>Visible</strong></p>\n`],
    [
      `---\ntitle: Keep me\n---\nAfter`,
      `<hr>\n<h2 id="title-keep-me">title: Keep me</h2>\n<p>After</p>\n`,
    ],
  ])(
    `preserves Markdown separators in descriptions: %s`,
    async (description, expected) => {
      const data = { example: { description, description_html: `` } }
      await render_data_markdown(data, `/repo/data/datasets.yml`)
      const rendered = document.createElement(`div`)
      rendered.innerHTML = data.example.description_html
      rendered
        .querySelectorAll(`[data-heading-anchor]`)
        .forEach((anchor) => anchor.remove())
      expect(rendered.innerHTML).toBe(expected)
    },
  )

  it.each([
    [`/repo/matbench_discovery/data-files.yml`, {}, `_links must be a string`],
    [`/repo/data/datasets.yml`, { bad: { description: 7 } }, `bad.description`],
    [
      `/repo/data/datasets.yml`,
      { bad: { description: ``, notes: { Access: 7 } } },
      `bad.notes.Access`,
    ],
    [`/repo/models/example/model.yml`, { notes: `bad` }, `notes: expected a mapping`],
  ])(`rejects invalid Markdown data in %s`, async (filename, data, message) => {
    await expect(render_data_markdown(data, filename)).rejects.toThrow(message)
  })
  it.each([
    [`Test String`, `test-string`],
    [`test_string`, `test-string`],
    [`Test__Multiple   Spaces`, `test-multiple-spaces`],
  ])(`assigns dataset slug '%s' → '%s'`, async (input, expected) => {
    const data = { [input]: { description: `` } }
    await render_data_markdown(data, `/repo/data/datasets.yml`)
    expect(data[input]).toEqual({ description: ``, description_html: ``, slug: expected })
  })
})

describe(`arr_to_str`, () => {
  it.each([
    [null, `n/a`],
    [undefined, `n/a`],
    [``, `n/a`],
    [[`a`, `b`, `c`], `a, b, c`],
    [123, `123`],
    [0, `0`],
    [true, `true`],
    [false, `false`],
  ] as const)(`converts %s → '%s'`, (input, expected) => {
    expect(arr_to_str(input as Parameters<typeof arr_to_str>[0])).toBe(expected)
  })
})

describe(`format_date`, () => {
  // covers string and numeric (timestamp) inputs; explicit times avoid timezone boundary issues
  it.each<[string | number, string, RegExp, RegExp]>([
    [`2023-05-15T12:00:00`, `2023`, /May/, /15/],
    [`2023-12-25T12:00:00`, `2023`, /Dec/, /25/],
    [new Date(`2023-12-25T12:00:00`).getTime(), `2023`, /Dec/, /25/],
  ])(`formats %s into a valid localized date`, (input, year, month, day) => {
    const result = format_date(input)
    expect(result).not.toBe(`Invalid Date`)
    expect(result).toContain(year)
    expect(result).toMatch(month)
    expect(result).toMatch(day)
  })
})

describe(`valid_query_param`, () => {
  it.each([
    [`F1`, `F1`],
    [`rmsd`, `rmsd`],
    [`model_params`, `model_params`],
    [`constructor`, `CPS`],
    [`__proto__`, `CPS`],
    [`md_time_multiplier`, `CPS`], // Table-only metric, not a scatter axis
    [``, `CPS`],
    [null, `CPS`],
  ])(`validates scatter axis %s against the option catalog`, (value, expected) => {
    const params = new URLSearchParams(value === null ? {} : { x: value })
    expect(valid_query_param(params, `x`, `CPS`, scatter_options_by_key)).toBe(expected)
  })
})

describe(`sync_url_params`, () => {
  it(`preserves unrelated params and omits defaults`, () => {
    history.replaceState(null, ``, `/benchmarks/md?keep=1&x=old&y=default#matrix`)

    sync_url_params(
      [
        [`x`, `new`, `default`],
        [`y`, `default`, `default`],
      ],
      {},
    )

    expect(location.search).toBe(`?keep=1&x=new`)
    expect(location.hash).toBe(`#matrix`)
  })

  it(`does not replace URL when params are unchanged`, () => {
    history.replaceState(null, ``, `/benchmarks/md?x=force_rmse`)
    const replace_spy = vi.spyOn(history, `replaceState`)

    sync_url_params([[`x`, `force_rmse`]], {})

    expect(replace_spy).not.toHaveBeenCalled()
  })

  it(`serializes writes so concurrent bindings keep each other's params`, async () => {
    // Kit's shallow goto() awaits route resolution before it writes history
    const deferred_goto = async (url: string | URL): Promise<void> => {
      await Promise.resolve()
      history.replaceState(null, ``, url)
    }
    vi.mocked(goto)
      .mockClear()
      .mockImplementationOnce(deferred_goto)
      .mockImplementationOnce(deferred_goto)

    sync_url_params([[`weights`, `0.5,0.5,0`]], {}) // root layout binding
    sync_url_params([[`train`, `MPtrj`]], {}) // page binding, same effect flush

    await vi.waitFor(() => expect(goto).toHaveBeenCalledTimes(2))
    await vi.waitFor(() => expect(location.search).toBe(`?weights=0.5,0.5,0&train=MPtrj`))
  })

  it(`keeps commas unencoded for human-readable weights params`, () => {
    history.replaceState(null, ``, `/`)

    sync_url_params([[`weights`, `0.579,0.35,0.071`]], {})

    expect(location.search).toBe(`?weights=0.579,0.35,0.071`)
    // round-trip: literal commas parse back to the same value
    expect(new URLSearchParams(location.search).get(`weights`)).toBe(`0.579,0.35,0.071`)
  })
})
