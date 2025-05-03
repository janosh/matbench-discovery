import {
  format_power_ten,
  format_property_path,
  get_format,
  get_org_logo,
} from '$lib/labels'
import { describe, expect, test } from 'vitest'

describe(`format_power_ten`, () => {
  test.each([
    {
      input: `1.23e-4`,
      expected: `1.23×10<sup>-4</sup>`,
      description: `negative exponent`,
    },
    {
      input: `5.67e+8`,
      expected: `5.67×10<sup>8</sup>`,
      description: `positive exponent with plus sign`,
    },
    {
      input: `9.01e12`,
      expected: `9.010<sup>12</sup>`,
      description: `positive exponent without plus sign`,
    },
    {
      input: `1e6`,
      expected: `10<sup>6</sup>`,
      description: `simplifies 1×10 to just 10`,
    },
    {
      input: `some text 1.23e-4 more text`,
      expected: `some text 1.23×10<sup>-4</sup> more text`,
      description: `works within text`,
    },
    {
      input: `1.23E-4`, // uppercase E
      expected: `1.23×10<sup>-4</sup>`,
      description: `uppercase E notation`,
    },
    {
      input: `multiple 1.2e3 and 4.5e-6 values`,
      expected: `multiple 1.2×10<sup>3</sup> and 4.5×10<sup>-6</sup> values`,
      description: `multiple scientific notations`,
    },
  ])(`formats '$input' to '$expected' ($description)`, ({ input, expected }) => {
    expect(format_power_ten(input)).toBe(expected)
  })

  test(`returns input string when no scientific notation is present`, () => {
    const input = `just a regular string with numbers 123.456`
    expect(format_power_ten(input)).toBe(input)
  })

  test(`handles edge cases correctly`, () => {
    // Empty string
    expect(format_power_ten(``)).toBe(``)

    // Edge case with 1×10
    expect(format_power_ten(`1×10<sup>3</sup>`)).toBe(`10<sup>3</sup>`)

    // Already formatted
    const formatted = `2.5×10<sup>-3</sup>`
    expect(format_power_ten(formatted)).toBe(formatted)
  })
})

describe(`format_property_path`, () => {
  test(`handles direct properties`, () => {
    expect(format_property_path(`model_params`)).toBe(`Number of model parameters`)
    expect(format_property_path(`date_added`)).toBe(`Date Added`)
    expect(format_property_path(`n_estimators`)).toBe(`Number of estimators`)
  })

  test(`returns the original path for unknown direct properties`, () => {
    expect(format_property_path(`unknown_property`)).toBe(`unknown property`)
  })

  test(`formats discovery metrics correctly`, () => {
    expect(format_property_path(`discovery.unique_prototypes.F1`)).toBe(
      `Discovery > Unique Prototypes > F1 Score`,
    )
    expect(format_property_path(`discovery.full_test_set.RMSE`)).toBe(
      `Discovery > Full Test Set > RMSE`,
    )
    expect(format_property_path(`discovery.most_stable_10k.Precision`)).toBe(
      `Discovery > 10k Most Stable > Precision`,
    )
  })

  test(`formats hyperparameter paths correctly`, () => {
    expect(format_property_path(`hyperparams.learning_rate`)).toBe(
      `hyperparams > Learning rate`,
    )
    expect(format_property_path(`hyperparams.graph_construction_radius`)).toBe(
      `hyperparams > Graph construction radius r<sub>cut</sub>`,
    )
    expect(format_property_path(`hyperparams.max_neighbors`)).toBe(
      `hyperparams > Max number of neighbors during graph construction`,
    )
    expect(format_property_path(`hyperparams.custom_param`)).toBe(
      `hyperparams > custom param`,
    )
  })

  test(`formats geo_opt metrics with symprec pattern correctly`, () => {
    expect(format_property_path(`geo_opt.symprec=1e-5.rmsd`)).toBe(
      `Geometry Optimization > RMSD`,
    )
  })

  test(`formats phonon metrics correctly`, () => {
    expect(format_property_path(`phonons.kappa_103.κ_SRME`)).toBe(
      `Phonons > κ<sub>SRME</sub>`,
    )
    expect(format_property_path(`phonons.other.property`)).toBe(
      `Phonons > other > property`,
    )
  })

  test(`handles generic dotted paths correctly`, () => {
    expect(format_property_path(`category.subcategory.property`)).toBe(
      `category > subcategory > property`,
    )
  })

  test(`handles special formatting in path parts`, () => {
    expect(format_property_path(`category.value_1e-5.property`)).toBe(
      `category > value 10<sup>-5</sup> > property`,
    )
  })

  test(`handles edge cases correctly`, () => {
    // Empty path
    expect(format_property_path(``)).toBe(``)

    // Path with empty segments
    expect(format_property_path(`..`)).toBe(``)

    // Path with single dot
    expect(format_property_path(`.`)).toBe(``)
  })
})

describe(`get_format`, () => {
  test(`returns '.1f' for empty arrays`, () => {
    expect(get_format([])).toBe(`.1f`)
  })

  test.each([
    {
      description: `date values (timestamps)`,
      values: [new Date(`2022-01-01`).getTime(), new Date(`2022-06-01`).getTime()],
      expected: `%b %y`,
    },
    {
      description: `large values (>10000)`,
      values: [15000, 20000, 25000],
      expected: `.1s`,
    },
    {
      description: `average > 1000`,
      values: [950, 1050, 1100],
      expected: `.1s`,
    },
    {
      description: `very small decimal values`,
      values: [0.0001, 0.0005, 0.0009],
      expected: `.5f`,
    },
    {
      description: `values with large range`,
      values: [1, 500, 1500],
      expected: `.2s`,
    },
    {
      description: `integer values`,
      values: [1, 2, 3, 4, 5],
      expected: `d`,
    },
    {
      description: `decimal values with moderate range`,
      values: [0.1, 0.5, 1.2, 2.3],
      expected: `.2f`,
    },
  ])(`returns '$expected' for $description`, ({ values, expected }) => {
    expect(get_format(values)).toBe(expected)
  })

  test(`identifies dates correctly near threshold boundaries`, () => {
    // Just after Jan 1, 2000
    const early_dates = [946_684_800_001, 946_684_900_000]
    expect(get_format(early_dates)).toBe(`%b %y`)

    // Just before Jan 1, 2050
    const late_dates = [2_524_607_999_999, 2_524_607_900_000]
    expect(get_format(late_dates)).toBe(`%b %y`)

    // Outside the valid date range - should not be treated as dates
    const before_range = [946_684_700_000, 946_684_750_000] // Before Jan 1, 2000
    const after_range = [2_524_608_100_000, 2_524_608_200_000] // After Jan 1, 2050

    expect(get_format(before_range)).not.toBe(`%b %y`)
    expect(get_format(after_range)).not.toBe(`%b %y`)
  })

  test(`handles mixed integer and decimal values correctly`, () => {
    const mixed_values = [1, 2, 3.5, 4.2]
    expect(get_format(mixed_values)).toBe(`.2f`)
  })

  test(`handles negative values correctly`, () => {
    const negative_values = [-10, -5, -3, -1]
    expect(get_format(negative_values)).toBe(`d`)

    const negative_mixed = [-100, -50.5, 25, 75.5]
    expect(get_format(negative_mixed)).toBe(`.2f`)

    const negative_small = [-0.001, -0.005, -0.009]
    expect(get_format(negative_small)).toBe(`.5f`)
  })

  test(`handles extreme values correctly`, () => {
    const extreme_values = [Number.MIN_SAFE_INTEGER, 0, Number.MAX_SAFE_INTEGER]
    // With extreme range, we expect scientific notation with 1 significant digit
    expect(get_format(extreme_values)).toBe(`.1s`)
  })

  test(`edge case: nearly-integer values are properly detected`, () => {
    // Values that are very close to integers but not exactly
    const nearly_integers = [1.0000001, 2.0000000001, 3]
    // JavaScript's floating point precision means these are treated as integers
    // The epsilon is likely too small for Math.abs(val - Math.round(val)) < 1e-6 check
    expect(get_format(nearly_integers)).toBe(`d`)

    // Values that are integers within floating point precision
    const float_integers = [1.0, 2.0, 3.0]
    // These should be treated as integers
    expect(get_format(float_integers)).toBe(`d`)
  })
})

describe(`get_org_logo`, () => {
  test.each([
    [`Google DeepMind`, { name: `Google DeepMind`, src: `/logos/deepmind.svg` }],
    [`FAIR at Meta`, { name: `FAIR at Meta`, id: `icon-logo-meta` }],
    [`Some unknown university`, undefined],
    [
      `Massachusetts Institute of Technology, USA`,
      { name: `Massachusetts Institute of Technology`, src: `/logos/mit.svg` },
    ],
    [`DeePMD`, { name: `DeePMD`, src: `/logos/deepmd.svg` }],
  ])(`returns correct logo data for '%s'`, (input, expected) => {
    expect(get_org_logo(input)).toEqual(expected)
  })

  test(`returns undefined for empty or undefined input`, () => {
    expect(get_org_logo(``)).toBeUndefined()
    // @ts-expect-error testing undefined input
    expect(get_org_logo(undefined)).toBeUndefined()
    // @ts-expect-error testing null input
    expect(get_org_logo(null)).toBeUndefined()
  })
})
