import { format_power_ten, format_property_path, get_org_logo } from '$lib/labels'
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
      `Hyperparams > Learning rate`,
    )
    expect(format_property_path(`hyperparams.graph_construction_radius`)).toBe(
      `Hyperparams > Graph construction radius r<sub>cut</sub>`,
    )
    expect(format_property_path(`hyperparams.max_neighbors`)).toBe(
      `Hyperparams > Max number of neighbors during graph construction`,
    )
    expect(format_property_path(`hyperparams.custom_param`)).toBe(
      `Hyperparams > custom param`,
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

describe(`get_org_logo`, () => {
  test.each([
    [`Google DeepMind`, { name: `Google DeepMind`, src: `/logos/deepmind.svg` }],
    [`FAIR at Meta`, { name: `FAIR at Meta`, id: `icon:Meta` }],
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
