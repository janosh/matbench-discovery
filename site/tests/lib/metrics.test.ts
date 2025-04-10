import { model_is_compliant } from '$lib'
import { calculate_cps, DEFAULT_CPS_CONFIG } from '$lib/combined_perf_score'
import {
  calc_cell_color,
  format_date,
  format_train_set,
  get_geo_opt_property,
  make_combined_filter,
  targets_tooltips,
} from '$lib/metrics'
import type { TargetType } from '$lib/model-schema'
import type { CombinedMetricConfig, ModelData } from '$lib/types'
import { describe, expect, it, test, vi } from 'vitest'

const model_is_compliant_mock = vi.hoisted(() => vi.fn())
const get_pred_file_urls_mock = vi.hoisted(() => vi.fn().mockReturnValue([]))

const mock_datasets = vi.hoisted(() => ({
  'MP 2022': {
    title: `Materials Project 2022`,
    url: `https://materialsproject.org`,
    n_structures: 100_000,
    n_materials: 50_000,
  },
  MPtrj: {
    title: `Materials Project Trajectories`,
    url: `https://materialsproject.org/trajectories`,
    n_structures: 200_000,
    n_materials: 75_000,
  },
  'Custom Set': {
    title: `Custom Dataset`,
    url: `https://example.com`,
    n_structures: 10_000,
  },
}))

// Mock the $lib module
vi.mock(`$lib`, () => ({
  model_is_compliant: model_is_compliant_mock,
  get_pred_file_urls: get_pred_file_urls_mock,
  DATASETS: mock_datasets,
}))

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

describe(`format_long_date`, () => {
  it(`formats date in long format`, () => {
    // Use a fixed date for testing
    const date = `2023-05-15`

    // Mock Date to return consistent results
    const original_date = Date
    const mock_date = class extends Date {
      constructor(date_str?: string | number | Date) {
        if (date_str) {
          super(date_str as string | number)
        } else {
          super(`2023-05-15T12:00:00Z`) // Mock current date
        }
      }

      toLocaleDateString(): string {
        return `Monday, May 15, 2023`
      }
    }

    // Override Date constructor
    Object.defineProperty(globalThis, `Date`, {
      value: mock_date,
      writable: true,
    })

    expect(format_date(date)).toBe(`Monday, May 15, 2023`)

    // Restore original Date
    Object.defineProperty(globalThis, `Date`, {
      value: original_date,
      writable: true,
    })
  })
})

describe(`format_train_set`, () => {
  it.each([
    {
      case: `single training set`,
      input: [`MP 2022`],
      expected_contains: [
        `<span title="`,
        `data-sort-value="50000"`,
        `<a href="https://materialsproject.org" target="_blank" rel="noopener noreferrer">MP 2022</a>`,
        `50k`,
        `50,000 materials in training set`,
      ],
    },
    {
      case: `multiple training sets`,
      input: [`MP 2022`, `MPtrj`],
      expected_contains: [
        `<a href="https://materialsproject.org" target="_blank" rel="noopener noreferrer">MP 2022</a>`,
        `<a href="https://materialsproject.org/trajectories" target="_blank" rel="noopener noreferrer">MPtrj</a>`,
        `MP 2022</a>+<a`,
        `data-sort-value="125000"`,
        `125k`,
        `&#013;• Materials Project 2022: 50,000 materials`,
        `• Materials Project Trajectories: 75,000 materials`,
      ],
    },
    {
      case: `training set with materials and structures`,
      input: [`MP 2022`],
      expected_contains: [
        `50k <small>(100k)</small>`,
        `50,000 materials in training set (100,000 structures`,
      ],
    },
    {
      case: `training set without n_materials`,
      input: [`Custom Set`],
      expected_contains: [
        `data-sort-value="10000"`,
        `10k`,
        `<a href="https://example.com" target="_blank" rel="noopener noreferrer">Custom Set</a>`,
      ],
    },
  ])(`formats $case correctly`, ({ input, expected_contains, not_contains }) => {
    const result = format_train_set(input)

    // Check for expected content
    for (const content of expected_contains) {
      expect(result).toContain(content)
    }

    // Check for content that should not be present
    if (not_contains) {
      for (const content of not_contains) {
        expect(result).not.toContain(content)
      }
    }
  })

  it(`handles missing training sets gracefully with warnings`, () => {
    // Mock console.warn
    const console_spy = vi.spyOn(console, `warn`).mockImplementation(() => {})

    const result = format_train_set([`MP 2022`, `NonExistent`])

    // Should warn about missing training set with exact message
    expect(console_spy).toHaveBeenCalledWith(
      `Training set NonExistent not found in DATASETS`,
    )

    // Should still format the existing training set correctly
    expect(result).toContain(
      `<a href="https://materialsproject.org" target="_blank" rel="noopener noreferrer">MP 2022</a>`,
    )
    expect(result).toContain(`50k`)

    // Should not include the missing dataset name anywhere
    expect(result).not.toContain(`NonExistent`)

    console_spy.mockRestore()
  })

  it(`shows both materials and structures when they differ`, () => {
    const result = format_train_set([`MP 2022`])

    // Verify both material and structure counts are shown
    expect(result).toContain(`50k <small>(100k)</small>`)

    // Check tooltip includes both counts with proper formatting
    expect(result).toContain(`50,000 materials in training set (100,000 structures`)
  })

  it(`formats training sets without n_materials correctly using n_structures`, () => {
    const result = format_train_set([`Custom Set`])

    // Should use n_structures as n_materials
    expect(result).toContain(`data-sort-value="10000"`)
    expect(result).toContain(`10k`)

    // Should include Custom Dataset in the result
    expect(result).toContain(
      `<a href="https://example.com" target="_blank" rel="noopener noreferrer">Custom Set</a>`,
    )
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

describe(`create_combined_filter`, () => {
  it.each([
    {
      case: `user filter returns false`,
      model_filter_returns: false,
      show_energy: true,
      show_noncomp: true,
      model: { targets: `E` },
      expected_result: false,
      should_check_compliance: false,
    },
    {
      case: `energy model with show_energy=false`,
      model_filter_returns: true,
      show_energy: false,
      show_noncomp: true,
      model: { targets: `E` },
      expected_result: false,
      should_check_compliance: false,
    },
    {
      case: `energy model with show_energy=true`,
      model_filter_returns: true,
      show_energy: true,
      show_noncomp: true,
      model: { targets: `E` },
      expected_result: true,
      should_check_compliance: true,
      is_compliant: true,
    },
    {
      case: `force model`,
      model_filter_returns: true,
      show_energy: false,
      show_noncomp: true,
      model: { targets: `EF_G` },
      expected_result: true,
      should_check_compliance: true,
      is_compliant: true,
    },
    {
      case: `non-compliant model with show_noncomp=false`,
      model_filter_returns: true,
      show_energy: true,
      show_noncomp: false,
      model: { targets: `EF_G` },
      expected_result: false,
      should_check_compliance: true,
      is_compliant: false,
    },
    {
      case: `non-compliant model with show_noncomp=true`,
      model_filter_returns: true,
      show_energy: true,
      show_noncomp: true,
      model: { targets: `EF_G` },
      expected_result: true,
      should_check_compliance: true,
      is_compliant: false,
    },
    {
      case: `all conditions pass`,
      model_filter_returns: true,
      show_energy: true,
      show_noncomp: true,
      model: { targets: `EF_G` },
      expected_result: true,
      should_check_compliance: true,
      is_compliant: true,
    },
  ])(
    `$case`,
    ({
      model_filter_returns,
      show_energy,
      show_noncomp,
      model,
      expected_result,
      should_check_compliance,
      is_compliant,
    }) => {
      const mock_model_filter = vi.fn().mockReturnValue(model_filter_returns)

      if (should_check_compliance && is_compliant !== undefined) {
        model_is_compliant_mock.mockReturnValue(is_compliant)
      }

      const filter = make_combined_filter(mock_model_filter, show_energy, show_noncomp)

      expect(filter(model as ModelData)).toBe(expected_result)
      expect(mock_model_filter).toHaveBeenCalledWith(model)

      if (should_check_compliance) {
        expect(model_is_compliant).toHaveBeenCalledWith(model)
      } else {
        expect(model_is_compliant).not.toHaveBeenCalled()
      }
    },
  )
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
