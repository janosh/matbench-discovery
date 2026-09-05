import {
  format_power_ten,
  format_property_path,
  format_relative_time,
  get_org_logo,
} from '$lib/labels'
import { Meta, Microsoft } from 'svelte-widgets/icons'
import { describe, expect, it } from 'vitest'

describe(`format_power_ten`, () => {
  it.each([
    [`negative exponent`, `1.23e-4`, `1.23×10<sup>-4</sup>`],
    [`positive exponent with plus sign`, `5.67e+8`, `5.67×10<sup>8</sup>`],
    [`mantissa ending in 1 is not collapsed`, `9.01e12`, `9.01×10<sup>12</sup>`],
    [`simplifies 1×10`, `1e6`, `10<sup>6</sup>`],
    [`collapses pre-formatted 1×10`, `1×10<sup>3</sup>`, `10<sup>3</sup>`],
    [`already formatted`, `2.5×10<sup>-3</sup>`, `2.5×10<sup>-3</sup>`],
    [`empty string`, ``, ``],
    [
      `no scientific notation`,
      `just a regular string with numbers 123.456`,
      `just a regular string with numbers 123.456`,
    ],
    [
      `within text`,
      `some text 1.23e-4 more text`,
      `some text 1.23×10<sup>-4</sup> more text`,
    ],
    [`uppercase E notation`, `1.23E-4`, `1.23×10<sup>-4</sup>`],
    [
      `multiple scientific notations`,
      `multiple 1.2e3 and 4.5e-6 values`,
      `multiple 1.2×10<sup>3</sup> and 4.5×10<sup>-6</sup> values`,
    ],
  ])(`formats %s`, (_description, input, expected) => {
    expect(format_power_ten(input)).toBe(expected)
  })
})

describe(`format_property_path`, () => {
  it.each([
    // Direct properties
    [`model_params`, `Params`],
    [`dates.benchmark_added`, `dates > Date Added`],
    [`n_estimators`, `Estimators`],
    [`unknown_property`, `unknown property`],
    // Discovery metrics
    [`discovery.unique_prototypes.F1`, `Discovery > Unique Prototypes > F1`],
    [`discovery.full_test_set.RMSE`, `Discovery > Full Test Set > RMSE`],
    [`discovery.unique_prototypes.Precision`, `Discovery > Unique Prototypes > Prec`],
    // Hyperparams
    [`hyperparams.training.learning_rate`, `Hyperparams > training > LR`],
    [
      `hyperparams.architecture.graph_construction_radius`,
      `Hyperparams > architecture > r<sub>cut</sub>`,
    ],
    [
      `hyperparams.architecture.max_neighbors`,
      `Hyperparams > architecture > Max neighbors`,
    ],
    [
      `hyperparams.upstream_config.custom_param`,
      `Hyperparams > upstream config > custom param`,
    ],
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
  ])(`formats '%s' → '%s'`, (input: string, expected: string) => {
    expect(format_property_path(input)).toBe(expected)
  })
})

describe(`get_org_logo`, () => {
  const src_logo = (name: string, src: string) => ({ name, src })

  it.each([
    [`Google DeepMind`, src_logo(`Google DeepMind`, `/logos/google-deepmind.svg`)],
    [`FAIR at Meta`, { name: `FAIR at Meta`, icon: Meta }],
    [
      `Microsoft Research AI for Science`,
      { name: `Microsoft Research`, icon: Microsoft },
    ],
    [`Some unknown university`, undefined],
    [
      `Massachusetts Institute of Technology, USA`,
      src_logo(
        `Massachusetts Institute of Technology`,
        `/logos/massachusetts-institute-of-technology.svg`,
      ),
    ],
    [
      `Department of Energy Science, Sungkyunkwan University`,
      src_logo(`Sungkyunkwan University`, `/logos/sungkyunkwan-university.svg`),
    ],
    [
      `Aberystwyth University, UK`,
      src_logo(`Aberystwyth University`, `/logos/aberystwyth-university.svg`),
    ],
    [
      `Independent Researcher in Data-driven Materials Discovery, M.Sc in Materials Science and Engineering, Christian-Albrechts-Universität zu Kiel`,
      src_logo(
        `Christian-Albrechts-Universität zu Kiel`,
        `/logos/christian-albrechts-university-kiel.svg`,
      ),
    ],
    [
      `Rutherford Appleton Laboratory, UK`,
      src_logo(
        `Rutherford Appleton Laboratory`,
        `/logos/science-and-technology-facilities-council.svg`,
      ),
    ],
    [`DeePMD`, src_logo(`DeePMD`, `/logos/deepmd.svg`)],
  ])(`returns correct logo data for '%s'`, (input: string, expected: unknown) => {
    expect(get_org_logo(input)).toStrictEqual(expected)
  })

  it(`returns undefined for empty or undefined input`, () => {
    expect(get_org_logo(``)).toBeUndefined()
    // @ts-expect-error testing undefined input
    expect(get_org_logo()).toBeUndefined()
    // @ts-expect-error testing null input
    expect(get_org_logo(null)).toBeUndefined()
  })
})

describe(`format_relative_time`, () => {
  const now = `2024-01-15T14:30:00Z`
  it.each([
    [`2024-01-15T14:29:30Z`, `1 minute ago`], // sub-minute clamps to 1
    [`2024-01-15T14:25:00Z`, `5 minutes ago`],
    [`2024-01-15T13:30:00Z`, `1 hour ago`],
    [`2024-01-15 09:30:00`, `5 hours ago`], // no timezone -> treated as UTC
    [`2024-01-14T14:30:00Z`, `1 day ago`],
    [`2024-01-01`, `14 days ago`],
    [`2024-02-01T08:05:00Z`, `2024-02-01 08:05:00 UTC`], // future -> absolute UTC
    [new Date(`2024-01-13T14:30:00Z`), `2 days ago`],
  ])(`formats %s relative to ${now} as '%s'`, (date, expected) => {
    expect(format_relative_time(date, now)).toBe(expected)
  })

  it(`returns N/A for missing or invalid dates`, () => {
    expect(format_relative_time(undefined, now)).toBe(`N/A`)
    expect(format_relative_time(`not a date`, now)).toBe(`N/A`)
    expect(format_relative_time(now, `garbage`)).toBe(`N/A`)
  })
})
