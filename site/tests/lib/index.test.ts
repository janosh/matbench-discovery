import { arr_to_str, calculate_days_ago, format_date, slugify } from '$lib'
import { describe, expect, it } from 'vitest'

describe(`calculate_days_ago`, () => {
  it(`returns empty string for empty input`, () => {
    expect(calculate_days_ago(``)).toBe(``)
  })

  it(`returns positive days for dates in the past`, () => {
    // We can test with a date that's definitely in the past
    const one_year_ago = new Date()
    one_year_ago.setFullYear(one_year_ago.getFullYear() - 1)
    const days_ago = calculate_days_ago(one_year_ago.toISOString())

    // Value should be approximately 365, but we just check it's positive
    const days_ago_num = parseInt(days_ago, 10)
    expect(days_ago_num).toBeGreaterThan(300)
    expect(days_ago_num).toBeLessThan(400)
  })

  it(`returns negative days for dates in the future`, () => {
    // We can test with a date that's definitely in the future
    const one_year_future = new Date()
    one_year_future.setFullYear(one_year_future.getFullYear() + 1)
    const days_ago = calculate_days_ago(one_year_future.toISOString())

    // Value should be approximately -365, but we just check it's negative
    const days_ago_num = parseInt(days_ago, 10)
    expect(days_ago_num).toBeLessThan(-300)
    expect(days_ago_num).toBeGreaterThan(-400)
  })
})

describe(`slugify`, () => {
  it(`converts spaces to hyphens and lowercases text`, () => {
    expect(slugify(`Test String`)).toBe(`test-string`)
  })

  it(`converts underscores to hyphens`, () => {
    expect(slugify(`test_string`)).toBe(`test-string`)
  })

  it(`handles multiple spaces and underscores`, () => {
    expect(slugify(`Test__Multiple   Spaces`)).toBe(`test-multiple-spaces`)
  })
})

describe(`arr_to_str`, () => {
  it(`returns n/a for null or undefined`, () => {
    expect(arr_to_str(null)).toBe(`n/a`)
    expect(arr_to_str(undefined)).toBe(`n/a`)
  })

  it(`joins arrays with commas`, () => {
    expect(arr_to_str([`a`, `b`, `c`])).toBe(`a, b, c`)
  })

  it(`converts non-array values to strings`, () => {
    expect(arr_to_str(123)).toBe(`123`)
    expect(arr_to_str(true)).toBe(`true`)
  })
})

describe(`format_date`, () => {
  it(`formats date with proper locale and options`, () => {
    const date = `2023-05-15`
    const result = format_date(date)

    expect(result).toContain(`2023`)
    expect(result).toMatch(/May/)
    expect(result).toMatch(/15/)
  })

  it(`handles different date formats`, () => {
    const iso_date = `2023-12-25`
    const timestamp = new Date(`2023-12-25`).getTime()

    const result1 = format_date(iso_date)
    const result2 = format_date(timestamp)

    // Both should produce valid date strings (not "Invalid Date")
    expect(result1).not.toBe(`Invalid Date`)
    expect(result2).not.toBe(`Invalid Date`)

    // Both should contain the year
    expect(result1).toContain(`2023`)
    expect(result2).toContain(`2023`)

    // Both should contain December or 12
    expect(result1).toMatch(/Dec/)
    expect(result2).toMatch(/Dec/)

    // Both should contain the day
    expect(result1).toMatch(/25/)
    expect(result2).toMatch(/25/)
  })
})
