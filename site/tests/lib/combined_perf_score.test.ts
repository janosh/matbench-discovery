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
function make_config(f1: number, rmsd: number, kappa: number): CpsConfig {
  return {
    F1: { ...DEFAULT_CPS_CONFIG.F1, weight: f1 },
    RMSD: { ...DEFAULT_CPS_CONFIG.RMSD, weight: rmsd },
    κ_SRME: { ...DEFAULT_CPS_CONFIG.κ_SRME, weight: kappa },
  }
}

describe(`calculate_cps`, () => {
  it(`returns weighted average when all metrics are perfect`, () => {
    const config = make_config(1 / 3, 1 / 3, 1 / 3)
    expect(calculate_cps(1.0, 0, 0, config)).toBeCloseTo(1.0, 5)
  })

  it(`returns 0 when all normalized metrics are worst`, () => {
    const config = make_config(1 / 3, 1 / 3, 1 / 3)
    expect(calculate_cps(0, RMSD_BASELINE, 2, config)).toBeCloseTo(0, 5)
  })

  it.each([
    { metric: `F1`, f1: undefined, rmsd: 0.05, kappa: 0.5 },
    { metric: `RMSD`, f1: 0.8, rmsd: undefined, kappa: 0.5 },
    { metric: `κ_SRME`, f1: 0.8, rmsd: 0.05, kappa: undefined },
    { metric: `F1 (NaN)`, f1: NaN, rmsd: 0.05, kappa: 0.5 },
    { metric: `RMSD (NaN)`, f1: 0.8, rmsd: NaN, kappa: 0.5 },
    { metric: `κ_SRME (NaN)`, f1: 0.8, rmsd: 0.05, kappa: NaN },
  ])(
    `returns null when $metric is missing/NaN with non-zero weight`,
    ({ f1, rmsd, kappa }) => {
      const config = make_config(0.5, 0.3, 0.2)
      expect(calculate_cps(f1, rmsd, kappa, config)).toBeNull()
    },
  )

  it.each([
    { metric: `F1`, weights: [0, 0.5, 0.5], f1: undefined, rmsd: 0, kappa: 0 },
    { metric: `RMSD`, weights: [0.5, 0, 0.5], f1: 1.0, rmsd: undefined, kappa: 0 },
    { metric: `κ_SRME`, weights: [0.5, 0.5, 0], f1: 1.0, rmsd: 0, kappa: undefined },
  ])(
    `ignores missing $metric when its weight is 0`,
    ({ weights, f1, rmsd, kappa }) => {
      const config = make_config(weights[0], weights[1], weights[2])
      expect(calculate_cps(f1, rmsd, kappa, config)).toBeCloseTo(1.0, 5)
    },
  )

  it(`returns 0 when all weights are 0`, () => {
    expect(calculate_cps(0.8, 0.05, 0.5, make_config(0, 0, 0))).toBe(0)
  })

  describe(`normalization`, () => {
    it.each([
      { name: `F1`, value: 0.5, weights: [1, 0, 0], expected: 0.5 },
      { name: `F1 perfect`, value: 1.0, weights: [1, 0, 0], expected: 1.0 },
      { name: `RMSD perfect`, value: 0, weights: [0, 1, 0], expected: 1.0 },
      { name: `RMSD worst`, value: RMSD_BASELINE, weights: [0, 1, 0], expected: 0.0 },
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
    const sum = DEFAULT_CPS_CONFIG.F1.weight +
      DEFAULT_CPS_CONFIG.κ_SRME.weight +
      DEFAULT_CPS_CONFIG.RMSD.weight
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
})
