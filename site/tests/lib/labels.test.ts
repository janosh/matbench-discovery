import { DEFAULT_CDS_CONFIG, DEFAULT_CMDS_CONFIG } from '#lib/combined-scores.svelte.js'
import {
  DIATOMICS_METRICS,
  GEO_OPT_SYMMETRY_METRICS,
  MD_METRICS,
  format_power_ten,
  format_property_path,
  format_relative_time,
  get_org_logo,
  plain_label,
} from '#lib/labels.js'
import { format_num } from 'matterviz/labels'
import { Meta, Microsoft } from 'svelte-widgets/icons'
import { describe, expect, it } from 'vite-plus/test'

it.each([
  [`E<sub>MAE</sub>`, `E_MAE`],
  [`R<sup>2</sup>`, `R^2`],
  [`t/t<sub>min</sub>`, `t/t_min`],
  [`Σ<sub>=</sub> 10<sup>-2</sup>`, `Σ_= 10^-2`],
  [`<span title="x">CPS</span>`, `CPS`],
  [`Params`, `Params`],
])(`plain_label(%s) = %s`, (label, expected) => {
  expect(plain_label(label)).toBe(expected)
})

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
    [`dates.benchmark_added`, `dates > Added`],
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
    // substring match returns the canonical org name, not the full affiliation
    [
      `Microsoft Research AI for Science`,
      { name: `Microsoft Research`, icon: Microsoft },
    ],
    [
      `Department of Energy Science, Sungkyunkwan University`,
      src_logo(`Sungkyunkwan University`, `/logos/sungkyunkwan-university.svg`),
    ],
    // case-insensitive match
    [`deepmd-kit team`, src_logo(`DeePMD`, `/logos/deepmd.svg`)],
    [`Some unknown university`, undefined],
    [``, undefined],
  ])(`returns correct logo data for '%s'`, (input: string, expected: unknown) => {
    expect(get_org_logo(input)).toStrictEqual(expected)
  })
})

it.each([
  [`2024-01-15T14:25:00Z`, `5m ago`],
  [`2024-01-01`, `2w 14h 30m ago`], // date-only strings parse as UTC midnight
  [`2022-11-20`, `1y 1mo 3w ago`],
  [`2024-02-01T08:05:00Z`, `2w 2d 17h from now`],
  [undefined, `N/A`],
  [`not a date`, `N/A`],
])(`format_relative_time(%s) relative to 2024-01-15T14:30Z is '%s'`, (date, expected) => {
  expect(format_relative_time(date, Date.parse(`2024-01-15T14:30:00Z`))).toBe(expected)
})

// descriptions quote the live default weights (rounded) and the CDS element roster size
const { vdos_error, adf_error, pressure_error, run_time_sec } = DEFAULT_CMDS_CONFIG
const { accuracy, geometry, physicality, speed } = DEFAULT_CDS_CONFIG
it.each([
  [
    MD_METRICS.md_combined_score,
    [
      [vdos_error, `ΔvDOS (%)`],
      [adf_error, `ΔADF (%)`],
      [pressure_error, `ΔP (%)`],
      [run_time_sec, `log-scaled, %;`],
    ],
  ],
  [
    DIATOMICS_METRICS.diatomics_combined_score,
    [
      [accuracy, `MAE, %)`],
      [geometry, `vs PBE, %)`],
      [physicality, `tortuosity, %)`],
      [speed, `log-scaled, %;`],
    ],
  ],
] as const)(`$label description quotes default weights`, (label, snippets) => {
  for (const [{ weight }, snippet] of snippets)
    expect(label.description).toContain(
      snippet.replace(`%`, `${Math.round(weight * 100)}%`),
    )
})
it(`CDS description quotes the scored element count`, () => {
  expect(DIATOMICS_METRICS.diatomics_combined_score.description).toContain(
    `fraction of 87 benchmark elements`,
  )
})

it.each(Object.values(MD_METRICS))(
  `$key label points at metrics.md and has correct direction`,
  (label) => {
    // t/t_min is computed per table view in assemble_row_data, not read from models
    const cost_keys = [`md_run_time_sec`, `md_max_gpu_mem_gb`, `md_max_rss_gb`]
    const expected_path =
      label.key === `md_time_multiplier`
        ? undefined
        : `metrics.md${cost_keys.includes(label.key) ? `.cost` : ``}`
    expect(label.path).toBe(expected_path)
    // combined_score (CMDS) is a score (higher=better); the rest are errors (lower)
    expect(label.better).toBe(label.key === `combined_score` ? `higher` : `lower`)
  },
)

// symmetry fractions show at most one decimal (YAMLs store 4, some more)
it.each([
  [0.74731771, `74.7%`],
  [0.0498, `5%`],
  [0.0009, `0.1%`],
])(`symmetry columns format %s as %s`, (value, expected) => {
  for (const { format } of Object.values(GEO_OPT_SYMMETRY_METRICS)) {
    expect(format_num(value, format)).toBe(expected)
  }
})
