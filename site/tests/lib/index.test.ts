import { arr_to_str, data_files, DATASETS, MODELS, format_date } from '$lib'
import { render_data_markdown } from '../../scripts/markdown-data'
import {
  apply_weights_param,
  sort_from_query,
  weights_to_param,
  sync_url_params,
} from '$lib/url-state.svelte'
import { valid_query_param } from 'svelte-widgets/url-params'
import { describe, expect, it, vi } from 'vitest'

describe(`$lib data includes rendered YAML Markdown`, () => {
  it(`DATASETS entries expose computed slug and description_html`, () => {
    const entries = Object.entries(DATASETS)
    expect(entries.length).toBeGreaterThanOrEqual(20) // datasets.yml entry count
    for (const [key, dataset] of entries) {
      for (const field of [`slug`, `description_html`] as const)
        expect(dataset[field]?.length, `${key} missing ${field}`).toBeGreaterThan(0)
    }
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
      expect(data.example.description_html).toBe(expected)
    },
  )

  it.each([
    [`/repo/matbench_discovery/data-files.yml`, {}, `_links must be a string`],
    [`/repo/data/datasets.yml`, { bad: { description: 7 } }, `bad.description`],
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
    [`direct key`, new URLSearchParams({ key: `valid` }), `valid`],
    [`prototype key`, new URLSearchParams({ key: `constructor` }), `fallback`],
    [`absent key`, new URLSearchParams(), `fallback`],
  ])(`returns expected value for %s`, (_case_name, params, expected) => {
    expect(valid_query_param(params, `key`, `fallback`, { valid: true })).toBe(expected)
  })
})

describe(`sort_from_query`, () => {
  const valid_cols = new Set([`F1`, `combined_score`])
  it.each([
    [`valid sort`, `?sort=rdf_error&dir=asc`, undefined, `rdf_error`, `asc`],
    [`invalid dir`, `?sort=F1&dir=sideways`, undefined, `F1`, `desc`],
    [`absent sort`, ``, undefined, `combined_score`, `desc`],
    [`column in valid_columns kept`, `?sort=F1`, valid_cols, `F1`, `desc`],
    [`unknown column falls back`, `?sort=nope`, valid_cols, `combined_score`, `desc`],
  ] as const)(`returns expected state for %s`, (_name, query, cols, column, dir) => {
    expect(
      sort_from_query(
        new URLSearchParams(query),
        { column: `combined_score`, dir: `desc` },
        cols,
      ),
    ).toStrictEqual({ column, dir })
  })
})

describe(`sync_url_params`, () => {
  it(`preserves unrelated params and omits defaults`, () => {
    history.replaceState(null, ``, `/tasks/md?keep=1&x=old&y=default#matrix`)

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
    history.replaceState(null, ``, `/tasks/md?x=force_rmse`)
    const replace_spy = vi.spyOn(history, `replaceState`)

    sync_url_params([[`x`, `force_rmse`]], {})

    expect(replace_spy).not.toHaveBeenCalled()
  })

  it(`keeps commas unencoded for human-readable weights params`, () => {
    history.replaceState(null, ``, `/`)

    sync_url_params([[`weights`, `0.579,0.35,0.071`]], {})

    expect(location.search).toBe(`?weights=0.579,0.35,0.071`)
    // round-trip: literal commas parse back to the same value
    expect(new URLSearchParams(location.search).get(`weights`)).toBe(`0.579,0.35,0.071`)
  })
})

describe(`weights_to_param / apply_weights_param`, () => {
  const make_config = (weights: number[]) => ({
    F1: { weight: weights[0] },
    κ_SRME: { weight: weights[1] },
    RMSD: { weight: weights[2] },
  })
  const defaults = make_config([0.5, 0.4, 0.1])

  it.each([
    [`defaults serialize to empty string`, [0.5, 0.4, 0.1], ``],
    [`custom weights serialize in key order`, [0.7, 0.2, 0.1], `0.7,0.2,0.1`],
    [`weights are rounded to 3 decimals`, [1 / 3, 1 / 3, 1 / 3], `0.333,0.333,0.333`],
  ] as const)(`%s`, (_name, weights, expected) => {
    expect(weights_to_param(make_config([...weights]), defaults)).toBe(expected)
  })

  it.each([
    [`valid param`, `0.7,0.2,0.1`, [0.7, 0.2, 0.1]],
    [`unnormalized weights get normalized`, `2,1,1`, [0.5, 0.25, 0.25]],
    // custom weights are shared module state; a weights-less URL must reset them,
    // and so must a malformed one - keeping stale in-session weights would show
    // wrong scores and launder them back into a valid-looking URL via the sync effect
    [`null param resets to defaults`, null, [0.5, 0.4, 0.1]],
    [`wrong count resets to defaults`, `0.5,0.5`, [0.5, 0.4, 0.1]],
    [`negative weight resets to defaults`, `-1,1,1`, [0.5, 0.4, 0.1]],
    [`non-numeric resets to defaults`, `a,b,c`, [0.5, 0.4, 0.1]],
    [`all-zero resets to defaults`, `0,0,0`, [0.5, 0.4, 0.1]],
    [`overflowing sum resets to defaults`, `1e308,1e308,1e308`, [0.5, 0.4, 0.1]],
  ] as const)(`%s`, (_name, param, expected) => {
    const config = make_config([0.2, 0.3, 0.5]) // start from non-default weights
    apply_weights_param(param, config, defaults)
    const weights = Object.values(config).map(({ weight }) => weight)
    for (const [idx, weight] of weights.entries()) {
      expect(weight).toBeCloseTo(expected[idx], 10)
    }
  })

  it(`round-trips through serialize -> parse`, () => {
    const config = make_config([0.62, 0.25, 0.13])
    const param = weights_to_param(config, defaults)
    const restored = make_config([0.5, 0.4, 0.1])
    apply_weights_param(param, restored, defaults)
    for (const key of Object.keys(config) as (keyof typeof config)[]) {
      expect(restored[key].weight).toBeCloseTo(config[key].weight, 3)
    }
  })
})
