// Tests for combined_perf_score.svelte.ts
import {
  calculate_cps,
  CPS_CONFIG,
  type CpsConfig,
  DEFAULT_CPS_CONFIG,
} from '$lib/combined_perf_score.svelte'
import { RMSD_BASELINE } from '$lib/labels'
import { afterEach, beforeEach, describe, expect, it } from 'vitest'

// Helper to create a test config with custom weights
const make_config = (f1: number, rmsd: number, kappa: number): CpsConfig => ({
  F1: { ...DEFAULT_CPS_CONFIG.F1, weight: f1 },
  RMSD: { ...DEFAULT_CPS_CONFIG.RMSD, weight: rmsd },
  κ_SRME: { ...DEFAULT_CPS_CONFIG.κ_SRME, weight: kappa },
})

describe(`calculate_cps`, () => {
  it.each([
    { name: `all perfect`, f1: 1.0, rmsd: 0, kappa: 0, expected: 1 },
    { name: `all worst`, f1: 0, rmsd: RMSD_BASELINE, kappa: 2, expected: 0 },
    { name: `mixed`, f1: 0.75, rmsd: RMSD_BASELINE / 2, kappa: 0.5, expected: 0.6667 },
  ])(
    `returns weighted average for $name metrics with equal weights`,
    ({ f1, rmsd, kappa, expected }) => {
      const config = make_config(1 / 3, 1 / 3, 1 / 3)
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
    { metric: `every metric`, f1: undefined, rmsd: undefined, kappa: undefined },
  ])(
    `returns null when $metric is missing/NaN with non-zero weight`,
    ({ f1, rmsd, kappa }) => {
      const config = make_config(0.5, 0.3, 0.2)
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
    // zero-weight metrics are ignored even when present with out-of-range values
    {
      case: `out-of-range values`,
      weights: [1, 0, 0],
      vals: [0.8, -0.01, 3],
      expected: 0.8,
    },
    // F1 passes through unclamped
    {
      case: `F1 above 1`,
      weights: [1, 0, 0],
      vals: [1.2, undefined, undefined],
      expected: 1.2,
    },
  ])(`ignores zero-weight metrics: $case`, ({ weights, vals, expected }) => {
    const config = make_config(weights[0], weights[1], weights[2])
    expect(calculate_cps(vals[0], vals[1], vals[2], config)).toBeCloseTo(expected, 5)
  })

  it(`returns 0 when all weights are 0`, () => {
    expect(calculate_cps(0.8, 0.05, 0.5, make_config(0, 0, 0))).toBe(0)
  })

  it(`normalizes weights that do not sum to 1`, () => {
    // weights sum to 2 -> effective (0.5, 0.25, 0.25); only the perfect F1 contributes
    const config = make_config(1.0, 0.5, 0.5)
    expect(calculate_cps(1.0, RMSD_BASELINE, 2, config)).toBeCloseTo(0.5, 4)
  })

  describe(`normalization`, () => {
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
      { name: `κ_SRME clamps high`, value: 3, weights: [0, 0, 1], expected: 0.0 },
    ])(`$name normalizes correctly`, ({ value, weights, expected }) => {
      const config = make_config(weights[0], weights[1], weights[2])
      const args: [number | undefined, number | undefined, number | undefined] = [
        weights[0] > 0 ? value : undefined,
        weights[1] > 0 ? value : undefined,
        weights[2] > 0 ? value : undefined,
      ]
      expect(calculate_cps(...args, config)).toBeCloseTo(expected, 5)
    })
  })

  it(`calculates correct CPS with default weights`, () => {
    // F1=0.8, RMSD=0.075 (midpoint), κ_SRME=0.5 -> 0.75
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
      (acc, { weight }) => acc + weight,
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
    // guards Reset + ?weights= URL persistence, both of which diff against defaults
    CPS_CONFIG.F1.weight = 0.9
    expect(DEFAULT_CPS_CONFIG.F1.weight).toBe(0.5)
  })
})
