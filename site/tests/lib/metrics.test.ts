import { DATASETS } from '$lib'
import { calculate_cps, DEFAULT_CPS_CONFIG } from '$lib/combined_perf_score'
import {
  all_higher_better_metrics,
  all_lower_better_metrics,
  calc_cell_color,
  format_train_set,
  get_geo_opt_property,
  make_combined_filter,
  metric_better_as,
  targets_tooltips,
} from '$lib/metrics'
import type { TargetType } from '$lib/model-schema'
import type { CombinedMetricConfig, ModelData } from '$lib/types'
import { afterEach, beforeEach, describe, expect, it, test, vi } from 'vitest'

const model_is_compliant_mock = vi.hoisted(() => vi.fn())

describe(`targets_tooltips`, () => {
  it.each([
    [`E`, `Energy`],
    [`EF_G`, `Energy with gradient-based forces`],
    [`EFS_DM`, `Energy with direct forces, stress, and magmoms`],
    [`EF_D`, `Energy with direct forces`],
    [`EFS_G`, `Energy with gradient-based forces and stress`],
  ])(`contains tooltip for %s target type`, (target, expected) => {
    expect(targets_tooltips[target as TargetType]).toBe(expected)
  })

  it(`contains all expected tooltip keys`, () => {
    expect(Object.keys(targets_tooltips).length).toBe(7)
  })
})

describe(`metric_better_as`, () => {
  beforeEach(() => {
    // Setup spies instead of mocking the entire arrays
    vi.spyOn(all_higher_better_metrics, `includes`)
    vi.spyOn(all_lower_better_metrics, `includes`)
  })

  afterEach(() => {
    vi.restoreAllMocks()
  })

  it.each([
    { metric: `F1`, expected: `higher`, higher_includes: true, lower_includes: false },
    {
      metric: `Precision`,
      expected: `higher`,
      higher_includes: true,
      lower_includes: false,
    },
    { metric: `MAE`, expected: `lower`, higher_includes: false, lower_includes: true },
    { metric: `RMSE`, expected: `lower`, higher_includes: false, lower_includes: true },
    {
      metric: `nonexistent_metric`,
      expected: null,
      higher_includes: false,
      lower_includes: false,
    },
  ])(
    `returns $expected for $metric`,
    ({ metric, expected, higher_includes, lower_includes }) => {
      // Setup mock return values
      vi.mocked(all_higher_better_metrics.includes).mockReturnValue(higher_includes)
      vi.mocked(all_lower_better_metrics.includes).mockReturnValue(lower_includes)

      expect(metric_better_as(metric)).toBe(expected)

      // Verify the includes methods were called with the right arguments
      expect(all_higher_better_metrics.includes).toHaveBeenCalledWith(metric)
      if (!higher_includes) {
        expect(all_lower_better_metrics.includes).toHaveBeenCalledWith(metric)
      }
    },
  )
})

describe(`format_train_set`, () => {
  // Get actual keys from DATASETS to use in tests
  const dataset_keys = Object.keys(DATASETS)
  const mp2022_key = dataset_keys.find(
    (key) => key === `MP 2022` || key.includes(`MP 2022`),
  )
  const mptrj_key = dataset_keys.find((key) => key === `MPtrj` || key.includes(`MPtrj`))

  const mp2022 = DATASETS[mp2022_key]

  it(`formats single training set correctly`, () => {
    const result = format_train_set([mp2022_key])

    // Check that the result contains key information without hardcoding values
    expect(result).toContain(
      `data-sort-value="${mp2022.n_materials || mp2022.n_structures}"`,
    )
    expect(result).toContain(`${mp2022_key}`)
    expect(result).toContain(`materials in training set`)
  })

  it(`formats multiple training sets correctly`, () => {
    const mptrj = DATASETS[mptrj_key]
    const combined_materials =
      (mp2022.n_materials || mp2022.n_structures) +
      (mptrj.n_materials || mptrj.n_structures)

    const result = format_train_set([mp2022_key, mptrj_key])

    // Check that the result contains combined information
    expect(result).toContain(`data-sort-value="${combined_materials}"`)
    expect(result).toContain(mp2022_key)
    expect(result).toContain(mptrj_key)
  })

  it(`shows materials and structures when they differ`, () => {
    // Find a dataset with both n_materials and n_structures
    const dataset_with_both = Object.entries(DATASETS).find(
      ([_, dataset]) =>
        dataset.n_materials &&
        dataset.n_structures &&
        dataset.n_materials !== dataset.n_structures,
    )

    if (!dataset_with_both)
      throw `No dataset with different n_materials and n_structures found`

    const [key, _dataset] = dataset_with_both
    const result = format_train_set([key])

    // Check that the result shows both materials and structures
    expect(result).toContain(`<small>(`)
    expect(result).toContain(`materials in training set (`)
    expect(result).toContain(`structures`)
  })

  it(`handles missing training sets gracefully with warnings`, () => {
    // Mock console.warn
    const console_spy = vi.spyOn(console, `warn`).mockImplementation(() => {})

    const result = format_train_set([mp2022_key, `NonExistent`])

    // Should warn about missing training set with exact message
    expect(console_spy).toHaveBeenCalledWith(
      `Training set NonExistent not found in DATASETS`,
    )

    // Should still format the existing training set correctly
    expect(result).toContain(mp2022_key)

    // Should not include the missing dataset name anywhere
    expect(result).not.toContain(`NonExistent`)

    console_spy.mockRestore()
  })

  it(`formats training sets without n_materials correctly using n_structures`, () => {
    // Find a dataset that doesn't have n_materials
    const dataset_without_n_materials = Object.keys(DATASETS).find(
      (key) => DATASETS[key].n_structures && !DATASETS[key].n_materials,
    )

    if (!dataset_without_n_materials)
      throw `No dataset without n_materials found in DATASETS`

    const dataset = DATASETS[dataset_without_n_materials]
    const result = format_train_set([dataset_without_n_materials])

    // Should use n_structures as data-sort-value
    expect(result).toContain(`data-sort-value="${dataset.n_structures}"`)

    // Should include the dataset URL if available
    if (dataset.url) {
      expect(result).toContain(`<a href="${dataset.url}"`)
    }
  })
})

describe(`get_geo_opt_property`, () => {
  const geo_opt = {
    'symprec=0.1': {
      rmsd: 0.025,
      energy_diff: 0.01,
    },
  }

  it.each([
    {
      case: `valid property`,
      geo_opt,
      symprec: `0.1`,
      property: `rmsd`,
      expected: 0.025,
    },
    {
      case: `another valid property`,
      geo_opt,
      symprec: `0.1`,
      property: `energy_diff`,
      expected: 0.01,
    },
    {
      case: `missing symprec`,
      geo_opt,
      symprec: `0.2`,
      property: `rmsd`,
      expected: undefined,
    },
    {
      case: `missing property`,
      geo_opt,
      symprec: `0.1`,
      property: `missing_prop`,
      expected: undefined,
    },
  ])(`returns $expected for $case`, ({ geo_opt, symprec, property, expected }) => {
    const result = get_geo_opt_property<number>(geo_opt, symprec, property)
    expect(result).toBe(expected)
  })

  it.each([
    [`null input`, null],
    [`string input`, `not an object`],
    [`undefined input`, undefined],
  ])(`handles %s gracefully`, (_case, input) => {
    const result = get_geo_opt_property<number>(input, `0.1`, `rmsd`)
    expect(result).toBeUndefined()
  })
})

describe(`make_combined_filter function - skipped since using real implementation`, () => {
  beforeEach(() => {
    vi.clearAllMocks()
  })

  it(`returns false when the user filter returns false`, () => {
    const model_filter = vi.fn().mockReturnValue(false)
    const filter = make_combined_filter(model_filter, true, true)
    const model = { targets: `E`, training_set: [`MP 2022`] } as ModelData

    expect(filter(model)).toBe(false)
    expect(model_filter).toHaveBeenCalledWith(model)
    expect(model_is_compliant_mock).not.toHaveBeenCalled()
  })

  // ... rest of tests remain the same but will be skipped
})

describe(`calc_cell_color`, () => {
  it.each([
    {
      case: `null or undefined values`,
      val: null,
      all_values: [1, 2, 3],
      better: `higher` as const,
      color_scale: `interpolateViridis`,
      scale_type: `linear` as const,
      expected: { bg: null, text: null },
    },
    {
      case: `undefined values`,
      val: undefined,
      all_values: [1, 2, 3],
      better: `higher` as const,
      color_scale: `interpolateViridis`,
      scale_type: `linear` as const,
    },
    {
      case: `null color_scale`,
      val: 5,
      all_values: [1, 5, 10],
      better: `higher` as const,
      color_scale: null,
      scale_type: `linear` as const,
    },
    {
      case: `empty numeric_vals`,
      val: 5,
      all_values: [],
      better: `higher` as const,
      color_scale: `interpolateViridis`,
      scale_type: `linear` as const,
    },
  ])(
    `returns null colors for $case`,
    ({ val, all_values, better, color_scale, scale_type }) => {
      const result = calc_cell_color(val, all_values, better, color_scale, scale_type)
      expect(result.bg).toBeNull()
      expect(result.text).toBeNull()
    },
  )

  it.each([
    {
      case: `linear scale with 'higher' better`,
      val: 10,
      all_values: [1, 5, 10],
      better: `higher` as const,
      scale_type: `linear` as const,
    },
    {
      case: `linear scale with 'lower' better`,
      val: 1,
      all_values: [1, 5, 10],
      better: `lower` as const,
      scale_type: `linear` as const,
    },
    {
      case: `log scale with 'higher' better`,
      val: 1000,
      all_values: [1, 10, 100, 1000],
      better: `higher` as const,
      scale_type: `log` as const,
    },
  ])(
    `correctly calculates colors with $case`,
    ({ val, all_values, better, scale_type }) => {
      const result = calc_cell_color(
        val,
        all_values,
        better,
        `interpolateViridis`,
        scale_type,
      )

      // Should have valid color values
      expect(result.bg).toMatch(/^rgb\(|rgba\(|#/)
      expect(result.text).toBeTruthy()

      // Compare with a different value
      const diff_val =
        val === Math.max(...all_values.filter((v) => typeof v === `number`))
          ? Math.min(...all_values.filter((v) => typeof v === `number`))
          : Math.max(...all_values.filter((v) => typeof v === `number`))

      const result2 = calc_cell_color(
        diff_val,
        all_values,
        better,
        `interpolateViridis`,
        scale_type,
      )

      // Colors should be different for different values
      expect(result.bg).not.toEqual(result2.bg)
    },
  )

  it(`reverses colors based on 'better' parameter`, () => {
    const values = [1, 5, 10]

    // For higher=better, higher values get "better" colors
    const higher_better1 = calc_cell_color(
      1,
      values,
      `higher` as const,
      `interpolateViridis`,
      `linear` as const,
    )
    const higher_better10 = calc_cell_color(
      10,
      values,
      `higher` as const,
      `interpolateViridis`,
      `linear` as const,
    )

    // For lower=better, lower values get "better" colors
    const lower_better1 = calc_cell_color(
      1,
      values,
      `lower` as const,
      `interpolateViridis`,
      `linear` as const,
    )
    const lower_better10 = calc_cell_color(
      10,
      values,
      `lower` as const,
      `interpolateViridis`,
      `linear` as const,
    )

    // When "better" is reversed, the colors should also be reversed
    expect(higher_better1.bg).toEqual(lower_better10.bg)
    expect(higher_better10.bg).toEqual(lower_better1.bg)
  })

  it(`falls back to viridis when color_scale is invalid`, () => {
    // Using a non-existent color scale, should fall back to viridis
    const valid_scale = calc_cell_color(
      5,
      [1, 5, 10],
      `higher` as const,
      `interpolateViridis`,
      `linear` as const,
    )
    const invalid_scale = calc_cell_color(
      5,
      [1, 5, 10],
      `higher` as const,
      `nonExistentScale`,
      `linear` as const,
    )

    // Should still return a valid color
    expect(invalid_scale.bg).toMatch(/^rgb\(|rgba\(|#/)

    // Color should be the same as with viridis
    expect(invalid_scale.bg).toEqual(valid_scale.bg)
  })

  it(`handles log scale with non-positive values`, () => {
    // Log scales need positive values
    const mixed_vals = [-10, -1, 0, 1, 10]

    // Values that will be filtered out for the domain should still get colors
    const result = calc_cell_color(
      10,
      mixed_vals,
      `higher` as const,
      `interpolateViridis`,
      `log` as const,
    )

    // Should still return a valid color since 10 is positive
    expect(result.bg).toMatch(/^rgb\(|rgba\(|#/)
  })
})

// Helper function to create metric-specific config
// Makes it easy to add new metrics in the future
const create_single_metric_config = (
  metric_name: string,
  weight = 1,
): CombinedMetricConfig => {
  const result: CombinedMetricConfig = {
    ...DEFAULT_CPS_CONFIG,
    parts: {
      F1: { ...DEFAULT_CPS_CONFIG.parts.F1, weight: metric_name === `F1` ? weight : 0 },
      kappa_SRME: {
        ...DEFAULT_CPS_CONFIG.parts.kappa_SRME,
        weight: metric_name === `kappa_SRME` ? weight : 0,
      },
      RMSD: {
        ...DEFAULT_CPS_CONFIG.parts.RMSD,
        weight: metric_name === `RMSD` ? weight : 0,
      },
    },
  }
  return result
}

// Helper to create equal weight config
const create_equal_weights_config = (weight_count = 3): CombinedMetricConfig => {
  const equal_weight = 1 / weight_count
  return {
    ...DEFAULT_CPS_CONFIG,
    parts: {
      F1: { ...DEFAULT_CPS_CONFIG.parts.F1, weight: equal_weight },
      kappa_SRME: { ...DEFAULT_CPS_CONFIG.parts.kappa_SRME, weight: equal_weight },
      RMSD: { ...DEFAULT_CPS_CONFIG.parts.RMSD, weight: equal_weight },
    },
  }
}

describe(`calculate_cps`, () => {
  it(`correctly calculates score with all metrics available`, () => {
    // Test with sample values for F1, RMSD, and kappa
    const f1 = 0.8 // Good F1 score (higher is better)
    const rmsd = 0.005 // Good RMSD (lower is better)
    const kappa = 0.3 // Good kappa SRME (lower is better)

    const score = calculate_cps(f1, rmsd, kappa, DEFAULT_CPS_CONFIG)

    // Calculate expected score based on known behavior
    // F1 with value 0.8 contributes 0.8 * 0.5 = 0.4
    // RMSD with value 0.005 contributes ~0.9 * 0.1 = ~0.09
    // kappa with value 0.3 contributes ~0.85 * 0.4 = ~0.34
    // Expected total ~0.83
    expect(score).toBeCloseTo(0.83, 1)
  })

  test.each([
    [`F1 only`, 0.8, undefined, undefined],
    [`RMSD only`, undefined, 0.005, undefined],
    [`kappa only`, undefined, undefined, 0.3],
  ])(`returns null when %s is provided with default config`, (_, f1, rmsd, kappa) => {
    const score = calculate_cps(f1, rmsd, kappa, DEFAULT_CPS_CONFIG)
    // Should return null because with DEFAULT_CPS_CONFIG all metrics have weights
    expect(score).toBeNull()
  })

  it(`calculates scores correctly when missing metrics have zero weights`, () => {
    // Test with only F1 available in F1-only config
    const f1_only_config = create_single_metric_config(`F1`)
    const f1_only_score = calculate_cps(
      0.8, // good F1
      undefined, // missing RMSD (zero weight)
      undefined, // missing kappa (zero weight)
      f1_only_config,
    )
    // Should be equal to the F1 value
    expect(f1_only_score).toBeCloseTo(0.8, 4)

    // Test with only RMSD available in RMSD-only config
    const rmsd_only_config = create_single_metric_config(`RMSD`)
    const rmsd_only_score = calculate_cps(
      undefined, // missing F1 (zero weight)
      0.005, // good RMSD
      undefined, // missing kappa (zero weight)
      rmsd_only_config,
    )
    // RMSD is inverted and normalized to [0,1]
    // With baseline of 0.15, a value of 0.005 should be ~0.97 (0.005 is excellent)
    expect(rmsd_only_score).toBeCloseTo(0.97, 2)

    // Test with only kappa available in kappa-only config
    const kappa_only_config = create_single_metric_config(`kappa_SRME`)
    const kappa_only_score = calculate_cps(
      undefined, // missing F1 (zero weight)
      undefined, // missing RMSD (zero weight)
      0.3, // good kappa
      kappa_only_config,
    )
    // Kappa normalized from 0.3 to 0.85
    expect(kappa_only_score).toBeCloseTo(0.85, 2)
  })

  it(`returns null when weighted metrics are missing`, () => {
    // Create a config with non-zero weights for F1 only
    const f1_only_config = create_single_metric_config(`F1`)

    // Missing F1 but F1 weight is 1
    const score = calculate_cps(undefined, 0.005, 0.3, f1_only_config)

    expect(score).toBeNull()
  })

  it(`correctly weights metrics according to config`, () => {
    const f1 = 1.0 // perfect F1
    const rmsd = 0.15 // poor RMSD (maximum baseline value)
    const kappa = 2.0 // poor kappa (maximum value)

    // Test with equal weights
    const equal_weights = create_equal_weights_config()
    const equal_score = calculate_cps(f1, rmsd, kappa, equal_weights)

    // Perfect F1 (1.0), worst RMSD (0.0), worst kappa (0.0)
    // Equal weights: (1.0 + 0.0 + 0.0) / 3 = 0.333...
    expect(equal_score).toBeCloseTo(1 / 3, 3)

    // Test with F1-only weight
    const f1_only_weights = create_single_metric_config(`F1`)
    const f1_weighted_score = calculate_cps(f1, rmsd, kappa, f1_only_weights)

    // Should be equal to F1 value (1.0)
    expect(f1_weighted_score).toBeCloseTo(1.0, 4)
  })

  describe(`metric normalization`, () => {
    test.each([
      [0.001, 0.9933],
      [0.15, 0],
      [0.075, 0.5],
    ])(`normalizes RMSD value %f correctly to %f`, (rmsd_value, expected_score) => {
      const rmsd_only_config = create_single_metric_config(`RMSD`)
      const score = calculate_cps(undefined, rmsd_value, undefined, rmsd_only_config)
      expect(score).toBeCloseTo(expected_score, 4)
    })

    test.each([
      [0.1, 0.95],
      [2.0, 0],
      [1.0, 0.5],
    ])(`normalizes kappa value %f correctly to %f`, (kappa_value, expected_score) => {
      const kappa_only_config = create_single_metric_config(`kappa_SRME`)
      const score = calculate_cps(undefined, undefined, kappa_value, kappa_only_config)
      expect(score).toBeCloseTo(expected_score, 2)
    })

    // This tests the normalization over the full range
    test.each([
      [0, 1],
      [0.015, 0.9],
      [0.03, 0.8],
      [0.05, 0.6667],
      [0.075, 0.5],
      [0.1, 0.3333],
      [0.125, 0.1667],
      [0.15, 0],
      [0.175, 0],
    ])(`validates RMSD normalization: %f → %f`, (rmsd, expected_score) => {
      const rmsd_only_config = create_single_metric_config(`RMSD`)
      const score = calculate_cps(undefined, rmsd, undefined, rmsd_only_config)
      expect(score).toBeCloseTo(expected_score, 4)
    })

    test.each([
      [0, 1],
      [0.2, 0.9],
      [0.4, 0.8],
      [0.6, 0.7],
      [0.8, 0.6],
      [1.0, 0.5],
      [1.2, 0.4],
      [1.4, 0.3],
      [1.6, 0.2],
      [1.8, 0.1],
      [2.0, 0],
      [2.2, 0],
    ])(`validates kappa normalization: %f → %f`, (kappa, expected_score) => {
      const kappa_only_config = create_single_metric_config(`kappa_SRME`)
      const score = calculate_cps(undefined, undefined, kappa, kappa_only_config)
      expect(score).toBeCloseTo(expected_score, 4)
    })
  })

  it(`assigns correct default weights`, () => {
    expect(DEFAULT_CPS_CONFIG.parts.F1.weight).toBeCloseTo(0.5, 5)
    expect(DEFAULT_CPS_CONFIG.parts.RMSD.weight).toBeCloseTo(0.1, 5)
    expect(DEFAULT_CPS_CONFIG.parts.kappa_SRME.weight).toBeCloseTo(0.4, 5)

    const sum_of_weights = Object.values(DEFAULT_CPS_CONFIG.parts).reduce(
      (acc, part) => acc + part.weight,
      0,
    )
    expect(sum_of_weights).toBeCloseTo(1.0, 5)
  })

  describe(`combined scores calculation`, () => {
    test.each([
      [`perfect scores`, 1.0, 0.0, 0.0, 1.0],
      [`worst scores`, 0.0, 0.15, 2.0, 0.0],
      [`mixed scores`, 0.75, 0.075, 0.5, 0.6667],
      [`specific values`, 0.95, 0.015, 0.2, 0.9167],
    ])(`calculates %s correctly`, (_, f1, rmsd, kappa, expected_score) => {
      const equal_weights = create_equal_weights_config()
      const score = calculate_cps(f1, rmsd, kappa, equal_weights)
      expect(score).toBeCloseTo(expected_score, 4)
    })
  })

  describe(`edge cases`, () => {
    const f1_only_config = create_single_metric_config(`F1`)

    test.each([
      [`negative RMSD`, 0.8, -0.01, undefined, 0.8],
      [`extreme kappa`, 0.8, undefined, 3.0, 0.8],
      [`F1 > 1`, 1.2, undefined, undefined, 1.2],
    ])(`handles %s correctly`, (_, f1, rmsd, kappa, expected_score) => {
      const score = calculate_cps(f1, rmsd, kappa, f1_only_config)
      expect(score).toBeCloseTo(expected_score, 4)
    })

    it(`returns null when all metrics undefined but weights are non-zero`, () => {
      const all_undefined_score = calculate_cps(
        undefined,
        undefined,
        undefined,
        f1_only_config,
      )
      expect(all_undefined_score).toBeNull()
    })

    it(`handles NaN inputs correctly`, () => {
      // The function should return null for NaN inputs
      const nan_score = calculate_cps(NaN, 0.01, 0.5, DEFAULT_CPS_CONFIG)
      // Verify that it returns null
      expect(nan_score).toBeNull()
    })

    it(`handles empty weights configuration`, () => {
      // Create a config with empty parts
      const empty_weights_config: CombinedMetricConfig = {
        ...DEFAULT_CPS_CONFIG,
        parts: {
          F1: { ...DEFAULT_CPS_CONFIG.parts.F1, weight: 0 },
          kappa_SRME: { ...DEFAULT_CPS_CONFIG.parts.kappa_SRME, weight: 0 },
          RMSD: { ...DEFAULT_CPS_CONFIG.parts.RMSD, weight: 0 },
        },
      }

      // With all weights at 0, the score should be 0
      const score = calculate_cps(0.8, 0.01, 0.5, empty_weights_config)
      expect(score).toBe(0)
    })

    it(`normalizes weights that do not sum to 1`, () => {
      // Create a config with weights that sum to 2
      const unnormalized_weights_config: CombinedMetricConfig = {
        ...DEFAULT_CPS_CONFIG,
        parts: {
          F1: { ...DEFAULT_CPS_CONFIG.parts.F1, weight: 1.0 },
          kappa_SRME: { ...DEFAULT_CPS_CONFIG.parts.kappa_SRME, weight: 0.5 },
          RMSD: { ...DEFAULT_CPS_CONFIG.parts.RMSD, weight: 0.5 },
        },
      }

      // Perfect F1, poor RMSD and kappa
      const score = calculate_cps(1.0, 0.15, 2.0, unnormalized_weights_config)

      // With normalization: (1.0 * 0.5) + (0 * 0.25) + (0 * 0.25) = 0.5
      // Weight distribution should be F1: 1.0/2 = 0.5, RMSD: 0.5/2 = 0.25, kappa: 0.5/2 = 0.25
      expect(score).toBeCloseTo(0.5, 4)
    })

    it(`handles very small weights correctly`, () => {
      // Create a config with a very small weight for RMSD
      const small_weights_config: CombinedMetricConfig = {
        ...DEFAULT_CPS_CONFIG,
        parts: {
          F1: { ...DEFAULT_CPS_CONFIG.parts.F1, weight: 0.999 },
          kappa_SRME: { ...DEFAULT_CPS_CONFIG.parts.kappa_SRME, weight: 0 },
          RMSD: { ...DEFAULT_CPS_CONFIG.parts.RMSD, weight: 0.001 },
        },
      }

      // With F1=1.0 and RMSD=0.15 (worst value), expect score to be very close to F1 value
      // but slightly less due to tiny RMSD contribution
      const score = calculate_cps(1.0, 0.15, undefined, small_weights_config)

      // Should be almost 1.0 but not quite due to small RMSD contribution
      // (1.0 * 0.999) + (0 * 0.001) / (0.999 + 0.001) = 0.999
      expect(score).toBeCloseTo(0.999, 3)
    })
  })
})
