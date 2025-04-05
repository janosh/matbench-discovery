import { calculate_combined_score, DEFAULT_CPS_CONFIG } from '$lib/metrics'
import type { CombinedMetricConfig } from '$lib/types'
import { describe, expect, it, test } from 'vitest'

describe(`Metrics`, () => {
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

  describe(`calculate_combined_score`, () => {
    it(`correctly calculates score with all metrics available`, () => {
      // Test with sample values for F1, RMSD, and kappa
      const f1 = 0.8 // Good F1 score (higher is better)
      const rmsd = 0.005 // Good RMSD (lower is better)
      const kappa = 0.3 // Good kappa SRME (lower is better)

      const score = calculate_combined_score(f1, rmsd, kappa, DEFAULT_CPS_CONFIG)

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
      const score = calculate_combined_score(f1, rmsd, kappa, DEFAULT_CPS_CONFIG)
      // Should return null because with DEFAULT_CPS_CONFIG all metrics have weights
      expect(score).toBeNull()
    })

    it(`calculates scores correctly when missing metrics have zero weights`, () => {
      // Test with only F1 available in F1-only config
      const f1_only_config = create_single_metric_config(`F1`)
      const f1_only_score = calculate_combined_score(
        0.8, // good F1
        undefined, // missing RMSD (zero weight)
        undefined, // missing kappa (zero weight)
        f1_only_config,
      )
      // Should be equal to the F1 value
      expect(f1_only_score).toBeCloseTo(0.8, 4)

      // Test with only RMSD available in RMSD-only config
      const rmsd_only_config = create_single_metric_config(`RMSD`)
      const rmsd_only_score = calculate_combined_score(
        undefined, // missing F1 (zero weight)
        0.005, // good RMSD
        undefined, // missing kappa (zero weight)
        rmsd_only_config,
      )
      // RMSD is inverted and normalized to [0,1]
      // With baseline of 0.03, a value of 0.005 should be ~0.83 (0.005 is excellent)
      expect(rmsd_only_score).toBeCloseTo(0.83, 2)

      // Test with only kappa available in kappa-only config
      const kappa_only_config = create_single_metric_config(`kappa_SRME`)
      const kappa_only_score = calculate_combined_score(
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
      const score = calculate_combined_score(undefined, 0.005, 0.3, f1_only_config)

      expect(score).toBeNull()
    })

    it(`correctly weights metrics according to config`, () => {
      const f1 = 1.0 // perfect F1
      const rmsd = 0.03 // poor RMSD (maximum baseline value)
      const kappa = 2.0 // poor kappa (maximum value)

      // Test with equal weights
      const equal_weights = create_equal_weights_config()
      const equal_score = calculate_combined_score(f1, rmsd, kappa, equal_weights)

      // Perfect F1 (1.0), worst RMSD (0.0), worst kappa (0.0)
      // Equal weights: (1.0 + 0.0 + 0.0) / 3 = 0.333...
      expect(equal_score).toBeCloseTo(1 / 3, 3)

      // Test with F1-only weight
      const f1_only_weights = create_single_metric_config(`F1`)
      const f1_weighted_score = calculate_combined_score(f1, rmsd, kappa, f1_only_weights)

      // Should be equal to F1 value (1.0)
      expect(f1_weighted_score).toBeCloseTo(1.0, 4)
    })

    describe(`metric normalization`, () => {
      test.each([
        [0.001, 0.9667],
        [0.03, 0],
        [0.015, 0.5],
      ])(`normalizes RMSD value %f correctly to %f`, (rmsd_value, expected_score) => {
        const rmsd_only_config = create_single_metric_config(`RMSD`)
        const score = calculate_combined_score(
          undefined,
          rmsd_value,
          undefined,
          rmsd_only_config,
        )
        expect(score).toBeCloseTo(expected_score, 4)
      })

      test.each([
        [0.1, 0.95],
        [2.0, 0],
        [1.0, 0.5],
      ])(`normalizes kappa value %f correctly to %f`, (kappa_value, expected_score) => {
        const kappa_only_config = create_single_metric_config(`kappa_SRME`)
        const score = calculate_combined_score(
          undefined,
          undefined,
          kappa_value,
          kappa_only_config,
        )
        expect(score).toBeCloseTo(expected_score, 2)
      })

      // This tests the normalization over the full range
      test.each([
        [0, 1],
        [0.003, 0.9],
        [0.006, 0.8],
        [0.01, 0.6667],
        [0.015, 0.5],
        [0.02, 0.3333],
        [0.025, 0.1667],
        [0.03, 0],
        [0.035, 0],
      ])(`validates RMSD normalization: %f → %f`, (rmsd, expected_score) => {
        const rmsd_only_config = create_single_metric_config(`RMSD`)
        const score = calculate_combined_score(
          undefined,
          rmsd,
          undefined,
          rmsd_only_config,
        )
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
        const score = calculate_combined_score(
          undefined,
          undefined,
          kappa,
          kappa_only_config,
        )
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
        [`worst scores`, 0.0, 0.03, 2.0, 0.0],
        [`mixed scores`, 0.75, 0.015, 0.5, 0.6667],
        [`specific values`, 0.95, 0.003, 0.2, 0.9167],
      ])(`calculates %s correctly`, (_, f1, rmsd, kappa, expected_score) => {
        const equal_weights = create_equal_weights_config()
        const score = calculate_combined_score(f1, rmsd, kappa, equal_weights)
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
        const score = calculate_combined_score(f1, rmsd, kappa, f1_only_config)
        expect(score).toBeCloseTo(expected_score, 4)
      })

      it(`returns null when all metrics undefined but weights are non-zero`, () => {
        const all_undefined_score = calculate_combined_score(
          undefined,
          undefined,
          undefined,
          f1_only_config,
        )
        expect(all_undefined_score).toBeNull()
      })

      it(`handles NaN inputs correctly`, () => {
        // The function should return null for NaN inputs
        const nan_score = calculate_combined_score(NaN, 0.01, 0.5, DEFAULT_CPS_CONFIG)
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
        const score = calculate_combined_score(0.8, 0.01, 0.5, empty_weights_config)
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
        const score = calculate_combined_score(
          1.0,
          0.03,
          2.0,
          unnormalized_weights_config,
        )

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

        // With F1=1.0 and RMSD=0.03 (worst value), expect score to be very close to F1 value
        // but slightly less due to tiny RMSD contribution
        const score = calculate_combined_score(1.0, 0.03, undefined, small_weights_config)

        // Should be almost 1.0 but not quite due to small RMSD contribution
        // (1.0 * 0.999) + (0 * 0.001) / (0.999 + 0.001) = 0.999
        expect(score).toBeCloseTo(0.999, 3)
      })
    })
  })
})
