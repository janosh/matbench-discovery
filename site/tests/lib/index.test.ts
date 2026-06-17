import { arr_to_str, data_files, DATASETS, format_date, slugify } from '$lib'
import { describe, expect, it } from 'vitest'

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
