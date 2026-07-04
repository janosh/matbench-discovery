import { arr_to_str, data_files, DATASETS, format_date, slugify } from '$lib'
import {
  sort_from_query,
  sync_url_params,
  valid_query_param,
} from '$lib/url-state.svelte'
import { describe, expect, it, vi } from 'vitest'

describe(`$lib data re-exports reflect index.ts mutations`, () => {
  it(`DATASETS entries expose computed slug and description_html`, () => {
    const entries = Object.entries(DATASETS)
    expect(entries.length).toBeGreaterThan(0)
    for (const [key, dataset] of entries) {
      expect(dataset.slug, `${key} missing slug`).toBe(slugify(key))
      expect(
        dataset.description_html?.length,
        `${key} missing description_html`,
      ).toBeGreaterThan(0)
    }
  })

  it(`data_files entries expose computed html`, () => {
    const entries = Object.entries(data_files).filter(
      ([key, entry]) => !key.startsWith(`_`) && typeof entry === `object`,
    )
    expect(entries.length).toBeGreaterThan(0)
    for (const [key, entry] of entries) {
      expect(
        (entry as { html?: string }).html?.length,
        `${key} missing html`,
      ).toBeGreaterThan(0)
    }
  })
})

describe(`slugify`, () => {
  it.each([
    [`Test String`, `test-string`],
    [`test_string`, `test-string`],
    [`Test__Multiple   Spaces`, `test-multiple-spaces`],
  ])(`converts '%s' → '%s'`, (input, expected) => {
    expect(slugify(input)).toBe(expected)
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
    history.replaceState(null, ``, `/tasks/md?keep=1&x=old&y=default`)

    sync_url_params(
      [
        [`x`, `new`, `default`],
        [`y`, `default`, `default`],
      ],
      {},
    )

    expect(location.search).toBe(`?keep=1&x=new`)
  })

  it(`does not replace URL when params are unchanged`, () => {
    history.replaceState(null, ``, `/tasks/md?x=force_rmse`)
    const replace_spy = vi.spyOn(history, `replaceState`)

    sync_url_params([[`x`, `force_rmse`]], {})

    expect(replace_spy).not.toHaveBeenCalled()
  })
})
