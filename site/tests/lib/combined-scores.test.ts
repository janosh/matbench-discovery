import type { ModelData } from '$lib'
import type {
  CpsConfig,
  CdsConfig,
  CdsValues,
  CmdsConfig,
} from '$lib/combined-scores.svelte'
import {
  calculate_cds,
  calculate_cps,
  CDS_COMPONENTS,
  CMDS_SPEED,
  CPS_CONFIG,
  DEFAULT_CDS_CONFIG,
  DEFAULT_CPS_CONFIG,
  update_models_cds,
  calculate_cmds,
  DEFAULT_CMDS_CONFIG,
  update_models_cmds,
} from '$lib/combined-scores.svelte'
import { RMSD_BASELINE } from '$lib/labels'
import { apply_weights_param, weights_to_param } from '$lib/url-state.svelte'
import { afterEach, beforeEach, describe, expect, it } from 'vitest'

const clone_config = (): CmdsConfig => structuredClone(DEFAULT_CMDS_CONFIG)
// geometric midpoint of the speed scale (log-scored), subscore exactly 0.5
const mid_speed = Math.sqrt(CMDS_SPEED.floor * CMDS_SPEED.baseline)
// config with speed zeroed and equal-weight observables = the pre-speed CMDS formula
const observables_only_config = (): CmdsConfig => {
  const config = clone_config()
  config.run_time_sec.weight = 0
  for (const key of [`adf_error`, `vdos_error`, `pressure_error`] as const) {
    config[key].weight = 1 / 3
  }
  return config
}
const clone_cds_config = (): CdsConfig => structuredClone(DEFAULT_CDS_CONFIG)
const make_cps_config = (f1: number, rmsd: number, kappa: number): CpsConfig => ({
  F1: { ...DEFAULT_CPS_CONFIG.F1, weight: f1 },
  RMSD: { ...DEFAULT_CPS_CONFIG.RMSD, weight: rmsd },
  κ_SRME: { ...DEFAULT_CPS_CONFIG.κ_SRME, weight: kappa },
})

// Every CDS component at the same fraction of its floor->baseline range (geometric
// interpolation for log-scaled components), so all subscores equal 1 - frac.
const cds_values_at = (frac: number): CdsValues =>
  Object.fromEntries(
    Object.values(CDS_COMPONENTS)
      .flat()
      .map(({ key, baseline, floor, log }) => [
        key,
        log ? floor * (baseline / floor) ** frac : floor + frac * (baseline - floor),
      ]),
  )
const perfect_cds_values = cds_values_at(0)
const halfway_cds_values = cds_values_at(0.5)
const baseline_cds_values = cds_values_at(1)

describe(`calculate_cps`, () => {
  it.each([
    { name: `all perfect`, f1: 1.0, rmsd: 0, kappa: 0, expected: 1 },
    { name: `all worst`, f1: 0, rmsd: RMSD_BASELINE, kappa: 2, expected: 0 },
    { name: `mixed`, f1: 0.75, rmsd: RMSD_BASELINE / 2, kappa: 0.5, expected: 0.6667 },
  ])(
    `returns weighted average for $name metrics with equal weights`,
    ({ f1, rmsd, kappa, expected }) => {
      const config = make_cps_config(1 / 3, 1 / 3, 1 / 3)
      expect(calculate_cps(f1, rmsd, kappa, config)).toBeCloseTo(expected, 4)
    },
  )

  it.each([
    { metric: `F1`, f1: undefined, rmsd: 0.05, kappa: 0.5 },
    { metric: `RMSD`, f1: 0.8, rmsd: undefined, kappa: 0.5 },
    { metric: `κ_SRME`, f1: 0.8, rmsd: 0.05, kappa: undefined },
    { metric: `F1 (NaN)`, f1: NaN, rmsd: 0.05, kappa: 0.5 },
    { metric: `RMSD (NaN)`, f1: 0.8, rmsd: NaN, kappa: 0.5 },
    { metric: `κ_SRME (NaN)`, f1: 0.8, rmsd: 0.05, kappa: NaN },
    { metric: `κ_SRME (Infinity)`, f1: 0.8, rmsd: 0.05, kappa: Infinity },
    { metric: `κ_SRME (negative)`, f1: 0.8, rmsd: 0.05, kappa: -1 },
    { metric: `κ_SRME (>2)`, f1: 0.8, rmsd: 0.05, kappa: 3 },
    { metric: `every metric`, f1: undefined, rmsd: undefined, kappa: undefined },
  ])(
    `returns null when $metric is missing/NaN with non-zero weight`,
    ({ f1, rmsd, kappa }) => {
      const config = make_cps_config(0.5, 0.3, 0.2)
      expect(calculate_cps(f1, rmsd, kappa, config)).toBeNull()
    },
  )

  it.each([
    { case: `missing F1`, weights: [0, 0.5, 0.5], vals: [undefined, 0, 0], expected: 1 },
    {
      case: `missing RMSD`,
      weights: [0.5, 0, 0.5],
      vals: [1, undefined, 0],
      expected: 1,
    },
    {
      case: `missing κ_SRME`,
      weights: [0.5, 0.5, 0],
      vals: [1, 0, undefined],
      expected: 1,
    },
    {
      case: `out-of-range values`,
      weights: [1, 0, 0],
      vals: [0.8, -0.01, 3],
      expected: 0.8,
    },
    {
      case: `F1 above 1`,
      weights: [1, 0, 0],
      vals: [1.2, undefined, undefined],
      expected: 1.2,
    },
  ])(`ignores zero-weight metrics: $case`, ({ weights, vals, expected }) => {
    const config = make_cps_config(weights[0], weights[1], weights[2])
    expect(calculate_cps(vals[0], vals[1], vals[2], config)).toBeCloseTo(expected, 5)
  })

  it(`returns null when all weights are 0 (score undefined, matches CMDS)`, () => {
    expect(calculate_cps(0.8, 0.05, 0.5, make_cps_config(0, 0, 0))).toBeNull()
  })

  it(`normalizes weights that do not sum to 1`, () => {
    const config = make_cps_config(1.0, 0.5, 0.5)
    expect(calculate_cps(1.0, RMSD_BASELINE, 2, config)).toBeCloseTo(0.5, 4)
  })

  it.each([
    { name: `F1`, value: 0.5, weights: [1, 0, 0], expected: 0.5 },
    { name: `F1 perfect`, value: 1.0, weights: [1, 0, 0], expected: 1.0 },
    { name: `RMSD perfect`, value: 0, weights: [0, 1, 0], expected: 1.0 },
    { name: `RMSD worst`, value: RMSD_BASELINE, weights: [0, 1, 0], expected: 0.0 },
    { name: `RMSD linear`, value: 0.03, weights: [0, 1, 0], expected: 0.8 },
    {
      name: `RMSD midpoint`,
      value: RMSD_BASELINE / 2,
      weights: [0, 1, 0],
      expected: 0.5,
    },
    {
      name: `RMSD clamps high`,
      value: RMSD_BASELINE + 0.1,
      weights: [0, 1, 0],
      expected: 0.0,
    },
    { name: `κ_SRME perfect`, value: 0, weights: [0, 0, 1], expected: 1.0 },
    { name: `κ_SRME worst`, value: 2, weights: [0, 0, 1], expected: 0.0 },
    { name: `κ_SRME linear`, value: 0.4, weights: [0, 0, 1], expected: 0.8 },
    { name: `κ_SRME midpoint`, value: 1, weights: [0, 0, 1], expected: 0.5 },
  ])(`$name normalizes correctly`, ({ value, weights, expected }) => {
    const config = make_cps_config(weights[0], weights[1], weights[2])
    const args: [number | undefined, number | undefined, number | undefined] = [
      weights[0] > 0 ? value : undefined,
      weights[1] > 0 ? value : undefined,
      weights[2] > 0 ? value : undefined,
    ]
    expect(calculate_cps(...args, config)).toBeCloseTo(expected, 5)
  })

  it(`calculates correct CPS with default weights`, () => {
    const result = calculate_cps(0.8, RMSD_BASELINE / 2, 0.5, DEFAULT_CPS_CONFIG)
    expect(result).toBeCloseTo(0.75, 5)
  })
})

describe(`DEFAULT_CPS_CONFIG`, () => {
  it(`has expected weights that sum to 1.0`, () => {
    expect(DEFAULT_CPS_CONFIG.F1.weight).toBe(0.5)
    expect(DEFAULT_CPS_CONFIG.κ_SRME.weight).toBe(0.4)
    expect(DEFAULT_CPS_CONFIG.RMSD.weight).toBe(0.1)
    const sum = Object.values(DEFAULT_CPS_CONFIG).reduce(
      (total, { weight }) => total + weight,
      0,
    )
    expect(sum).toBeCloseTo(1.0, 5)
  })
})

describe(`CPS_CONFIG reactivity`, () => {
  let original_weights: [number, number, number]

  beforeEach(() => {
    original_weights = [
      CPS_CONFIG.F1.weight,
      CPS_CONFIG.RMSD.weight,
      CPS_CONFIG.κ_SRME.weight,
    ]
  })

  afterEach(() => {
    ;[CPS_CONFIG.F1.weight, CPS_CONFIG.RMSD.weight, CPS_CONFIG.κ_SRME.weight] =
      original_weights
  })

  it(`can be modified and affects calculate_cps`, () => {
    CPS_CONFIG.F1.weight = 1.0
    CPS_CONFIG.RMSD.weight = 0
    CPS_CONFIG.κ_SRME.weight = 0
    expect(calculate_cps(0.6, 0, 0, CPS_CONFIG)).toBeCloseTo(0.6, 5)
  })

  it(`weight edits must not corrupt DEFAULT_CPS_CONFIG`, () => {
    CPS_CONFIG.F1.weight = 0.9
    expect(DEFAULT_CPS_CONFIG.F1.weight).toBe(0.5)
  })

  it(`round-trips custom weights through the ?weights= URL param`, () => {
    expect(weights_to_param(CPS_CONFIG, DEFAULT_CPS_CONFIG)).toBe(``)

    CPS_CONFIG.F1.weight = 0.7
    CPS_CONFIG.κ_SRME.weight = 0.2
    CPS_CONFIG.RMSD.weight = 0.1
    expect(weights_to_param(CPS_CONFIG, DEFAULT_CPS_CONFIG)).toBe(`0.7,0.2,0.1`)

    apply_weights_param(null, CPS_CONFIG, DEFAULT_CPS_CONFIG)
    expect(CPS_CONFIG.F1.weight).toBe(0.5)
    expect(weights_to_param(CPS_CONFIG, DEFAULT_CPS_CONFIG)).toBe(``)
  })
})

describe(`calculate_cmds`, () => {
  const perfect = { adf_error: 0, vdos_error: 0, pressure_error: 0 }
  it.each([
    [{ ...perfect, run_time_sec: CMDS_SPEED.floor }, 1],
    [
      {
        adf_error: 100,
        vdos_error: 100,
        pressure_error: 100,
        run_time_sec: CMDS_SPEED.baseline,
      },
      0,
    ],
    // default weights: 30% vdos, 20% adf, 20% speed, 30% pressure
    [
      { adf_error: 10, vdos_error: 20, pressure_error: 30, run_time_sec: mid_speed },
      0.3 * 0.8 + 0.2 * 0.9 + 0.2 * 0.5 + 0.3 * 0.7,
    ],
    // per-component subscores clamp to [0,1]: the >100% ADF error clamps to subscore
    // 0 and the faster-than-floor speed clamps to 1, so neither over/undershoots
    [{ adf_error: 150, vdos_error: 0, pressure_error: 0, run_time_sec: 1 }, 0.8],
  ])(`default-weight score for %o = %f`, (values, expected) => {
    expect(calculate_cmds(values, clone_config())).toBeCloseTo(expected, 10)
  })

  it(`zero-weight speed + equal observables reproduces the pre-speed formula`, () => {
    const errors = { adf_error: 10, vdos_error: 20, pressure_error: 30 }
    expect(calculate_cmds(errors, observables_only_config())).toBeCloseTo(0.8, 10)
  })

  it.each([
    [{ adf_error: 10, vdos_error: 20, run_time_sec: 100 }], // missing component
    [{ adf_error: 10, vdos_error: 20, pressure_error: NaN, run_time_sec: 100 }],
    [{ ...perfect, run_time_sec: 0 }], // log speed scale needs positive wall time
    [{ ...perfect }], // missing run_time_sec with non-zero speed weight
  ])(`returns null when a non-zero-weight component is invalid: %o`, (values) => {
    expect(calculate_cmds(values, clone_config())).toBeNull()
  })

  it(`skips missing components with zero weight`, () => {
    const config = observables_only_config()
    config.pressure_error.weight = 0
    // mean of (1 - 10/100) and (1 - 20/100) weighted 1/3 each = 0.85
    expect(calculate_cmds({ adf_error: 10, vdos_error: 20 }, config)).toBeCloseTo(0.85)
  })

  it(`returns null when all weights are zero (score undefined, not worst)`, () => {
    const config = clone_config()
    for (const key of Object.keys(config) as (keyof CmdsConfig)[]) {
      config[key].weight = 0
    }
    expect(calculate_cmds({}, config)).toBeNull()
  })
})

describe(`update_models_cmds`, () => {
  it(`writes CMDS into metrics.md and deletes it when components are missing`, () => {
    const models = [
      {
        metrics: {
          md: {
            adf_error: 10,
            vdos_error: 20,
            pressure_error: 30,
            run_time_sec: mid_speed,
          },
        },
      },
      // no run_time_sec: stale score must be removed under default (non-zero) speed weight
      {
        metrics: {
          md: {
            adf_error: 10,
            vdos_error: 20,
            pressure_error: 30,
            combined_score: 0.99,
          },
        },
      },
      { metrics: {} },
    ] as unknown as ModelData[]

    update_models_cmds(models, clone_config())

    const scored_md = models[0].metrics?.md as { combined_score?: number }
    expect(scored_md.combined_score).toBeCloseTo(
      0.3 * 0.8 + 0.2 * 0.9 + 0.2 * 0.5 + 0.3 * 0.7,
      10,
    )
    // incomplete components: stale score removed, not overwritten with NaN (a NaN
    // field would count as present in scatter model counts and export as 'NaN')
    expect(models[1].metrics?.md).not.toHaveProperty(`combined_score`)
    expect(models[2].metrics?.md).toBeUndefined()
  })
})

describe(`calculate_cds`, () => {
  it.each([
    [`perfect`, perfect_cds_values, 1],
    [`all at baseline`, baseline_cds_values, 0],
    [`all halfway`, halfway_cds_values, 0.5],
    [`beyond baseline`, { ...baseline_cds_values, pbe_energy_mae: 100 }, 0],
    [`below floor`, { ...perfect_cds_values, force_flips: 0.5 }, 1],
    [
      `pillar weighting`,
      { ...baseline_cds_values, pbe_energy_mae: 0, pbe_force_mae: 0 },
      4 / 9,
    ],
  ])(`default-weight score: %s = %f`, (_desc, values, expected) => {
    expect(calculate_cds(values, clone_cds_config())).toBeCloseTo(expected, 10)
  })

  it(`penalizes excluded benchmark elements through coverage`, () => {
    const values = {
      ...perfect_cds_values,
      excluded_formula_reasons: { [`He-He`]: `non-finite repulsive wall` },
    }
    expect(calculate_cds(values, clone_cds_config())).toBeCloseTo(86 / 87, 10)
  })

  it.each([
    [{ ...halfway_cds_values, tortuosity: undefined }],
    [{ ...halfway_cds_values, pbe_force_mae: Number.NaN }],
    [{ ...halfway_cds_values, pbe_force_mae: Number.POSITIVE_INFINITY }],
    [{ ...halfway_cds_values, run_time_sec: 0 }],
  ])(
    `returns null when a non-zero-weight pillar has an invalid component: %o`,
    (values) => {
      expect(calculate_cds(values as CdsValues, clone_cds_config())).toBeNull()
    },
  )

  it(`skips missing components in zero-weight pillars`, () => {
    const config = clone_cds_config()
    config.physicality.weight = 0
    const {
      energy_jump: _energy_jump,
      force_flips: _force_flips,
      tortuosity: _tortuosity,
      ...values
    } = halfway_cds_values

    expect(calculate_cds(values, config)).toBeCloseTo(0.5, 10)
  })

  it(`returns null when all weights are zero (score undefined, not worst)`, () => {
    const config = clone_cds_config()
    for (const key of Object.keys(config) as (keyof CdsConfig)[]) {
      config[key].weight = 0
    }
    expect(calculate_cds(halfway_cds_values, config)).toBeNull()
  })
})

describe(`update_models_cds`, () => {
  it(`writes CDS into metrics.diatomics and deletes it when components are missing`, () => {
    const models = [
      { metrics: { diatomics: { ...halfway_cds_values } } },
      { metrics: { diatomics: { pbe_energy_mae: 1, combined_score: 0.99 } } },
      { metrics: {} },
    ] as unknown as ModelData[]

    update_models_cds(models, clone_cds_config())

    const scored = models[0].metrics?.diatomics as { combined_score?: number }
    expect(scored.combined_score).toBeCloseTo(0.5, 10)
    expect(models[1].metrics?.diatomics).not.toHaveProperty(`combined_score`)
    expect(models[2].metrics?.diatomics).toBeUndefined()
  })
})
