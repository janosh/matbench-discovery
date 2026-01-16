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
  test.each([
    // Direct properties
    [`model_params`, `Params`],
    [`date_added`, `date added`],
    [`n_estimators`, `Estimators`],
    [`unknown_property`, `unknown property`],
    // Discovery metrics
    [`discovery.unique_prototypes.F1`, `Discovery > Unique Prototypes > F1`],
    [`discovery.full_test_set.RMSE`, `Discovery > Full Test Set > RMSE`],
    [`discovery.most_stable_10k.Precision`, `Discovery > 10k Most Stable > Prec`],
    // Hyperparameters
    [`hyperparams.learning_rate`, `Hyperparams > LR`],
    [`hyperparams.graph_construction_radius`, `Hyperparams > r<sub>cut</sub>`],
    [`hyperparams.max_neighbors`, `Hyperparams > Max neighbors`],
    [`hyperparams.custom_param`, `Hyperparams > custom param`],
    // Geo-opt metrics
    [`geo_opt.symprec=1e-5.rmsd`, `Geometry Optimization > RMSD`],
    // Phonon metrics
    [`phonons.kappa_103.κ_SRME`, `Phonons > κ<sub>SRME</sub>`],
    [`phonons.other.property`, `Phonons > other > property`],
    // Generic paths
    [`category.subcategory.property`, `category > subcategory > property`],
    [`category.value_1e-5.property`, `category > value 10<sup>-5</sup> > property`],
    // Edge cases
    [``, ``],
    [`..`, ``],
    [`.`, ``],
  ])(`formats '%s' → '%s'`, (input, expected) => {
    expect(format_property_path(input)).toBe(expected)
  })
})

describe(`get_org_logo`, () => {
  test.each([
    [`Google DeepMind`, { name: `Google DeepMind`, src: `/logos/deepmind.svg` }],
    [`FAIR at Meta`, { name: `FAIR at Meta`, id: `icon:LogoMeta` }],
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
