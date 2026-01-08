import { arr_to_str, calculate_days_ago, format_date, slugify } from '$lib'
import { describe, expect, it } from 'vitest'

describe(`calculate_days_ago`, () => {
  it(`returns empty string for empty input`, () => {
    expect(calculate_days_ago(``)).toBe(``)
  })

  it.each([
    { offset_years: -1, min: 300, max: 400, desc: `past date` },
    { offset_years: 1, min: -400, max: -300, desc: `future date` },
  ])(`handles $desc correctly`, ({ offset_years, min, max }) => {
    const date = new Date()
    date.setFullYear(date.getFullYear() + offset_years)
    const days_ago = parseInt(calculate_days_ago(date.toISOString()), 10)
    expect(days_ago).toBeGreaterThan(min)
    expect(days_ago).toBeLessThan(max)
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
  it.each(
    [
      [null, `n/a`],
      [undefined, `n/a`],
      [[`a`, `b`, `c`], `a, b, c`],
      [123, `123`],
      [true, `true`],
    ] as const,
  )(`converts %s → '%s'`, (input, expected) => {
    // @ts-expect-error testing various input types
    expect(arr_to_str(input)).toBe(expected)
  })
})

describe(`format_date`, () => {
  it(`formats date with proper locale and options`, () => {
    // Use a date with explicit time to avoid timezone issues
    const date = `2023-05-15T12:00:00`
    const result = format_date(date)

    expect(result).toContain(`2023`)
    expect(result).toMatch(/May/)
    expect(result).toMatch(/15/)
  })

  it(`handles different date formats`, () => {
    // Use explicit times to avoid timezone boundary issues
    const iso_date = `2023-12-25T12:00:00`
    const timestamp = new Date(`2023-12-25T12:00:00`).getTime()

    const result1 = format_date(iso_date)
    const result2 = format_date(timestamp)

    // Both should produce valid date strings (not "Invalid Date")
    expect(result1).not.toBe(`Invalid Date`)
    expect(result2).not.toBe(`Invalid Date`)

    // Both should contain the year
    expect(result1).toContain(`2023`)
    expect(result2).toContain(`2023`)

    // Both should contain December
    expect(result1).toMatch(/Dec/)
    expect(result2).toMatch(/Dec/)

    // Both should contain the day
    expect(result1).toMatch(/25/)
    expect(result2).toMatch(/25/)
  })
})
