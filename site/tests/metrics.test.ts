import {
  calculate_combined_score,
  DEFAULT_COMBINED_METRIC_CONFIG,
  F1_DEFAULT_WEIGHT,
  KAPPA_DEFAULT_WEIGHT,
  RMSD_DEFAULT_WEIGHT,
} from '$lib/metrics'
import { describe, expect, it } from 'vitest'

describe(`Metrics`, () => {
  describe(`calculate_combined_score`, () => {
    it(`correctly calculates score with all metrics available`, () => {
      // Test with sample values for F1, RMSD, and kappa
      const f1 = 0.8 // Good F1 score (higher is better)
      const rmsd = 0.05 // Good RMSD (lower is better)
      const kappa = 0.3 // Good kappa SRME (lower is better)

      const score = calculate_combined_score(
        f1,
        rmsd,
        kappa,
        DEFAULT_COMBINED_METRIC_CONFIG,
      )

      // Score should be a number between 0 and 1
      expect(score).toBeGreaterThan(0)
      expect(score).toBeLessThan(1)

      // For these good values, we expect a high score
      expect(score).toBeGreaterThan(0.7)
    })

    it(`returns NaN when metrics with non-zero weights are missing`, () => {
      // With DEFAULT_COMBINED_METRIC_CONFIG, all metrics have non-zero weights
      // So if any are missing, we should get NaN

      // Test with only F1 available
      const f1_only_score = calculate_combined_score(
        0.8, // good F1
        undefined, // missing RMSD
        undefined, // missing kappa
        DEFAULT_COMBINED_METRIC_CONFIG,
      )

      // Should return NaN because RMSD and kappa are missing but have weights
      expect(isNaN(f1_only_score)).toBe(true)

      // Test with only RMSD available
      const rmsd_only_score = calculate_combined_score(
        undefined, // missing F1
        0.05, // good RMSD
        undefined, // missing kappa
        DEFAULT_COMBINED_METRIC_CONFIG,
      )

      // Should return NaN
      expect(isNaN(rmsd_only_score)).toBe(true)

      // Test with only kappa available
      const kappa_only_score = calculate_combined_score(
        undefined, // missing F1
        undefined, // missing RMSD
        0.3, // good kappa
        DEFAULT_COMBINED_METRIC_CONFIG,
      )

      // Should return NaN
      expect(isNaN(kappa_only_score)).toBe(true)
    })

    it(`calculates scores correctly when missing metrics have zero weights`, () => {
      // Create configs where only one metric has weight
      const f1_only_config = {
        ...DEFAULT_COMBINED_METRIC_CONFIG,
        weights: [
          { metric: `F1`, value: 1, display: `F1`, description: `` },
          { metric: `RMSD`, value: 0, display: `RMSD`, description: `` },
          { metric: `kappa_SRME`, value: 0, display: `kappa`, description: `` },
        ],
      }

      const rmsd_only_config = {
        ...DEFAULT_COMBINED_METRIC_CONFIG,
        weights: [
          { metric: `F1`, value: 0, display: `F1`, description: `` },
          { metric: `RMSD`, value: 1, display: `RMSD`, description: `` },
          { metric: `kappa_SRME`, value: 0, display: `kappa`, description: `` },
        ],
      }

      const kappa_only_config = {
        ...DEFAULT_COMBINED_METRIC_CONFIG,
        weights: [
          { metric: `F1`, value: 0, display: `F1`, description: `` },
          { metric: `RMSD`, value: 0, display: `RMSD`, description: `` },
          { metric: `kappa_SRME`, value: 1, display: `kappa`, description: `` },
        ],
      }

      // Test with only F1 available in F1-only config
      const f1_only_score = calculate_combined_score(
        0.8, // good F1
        undefined, // missing RMSD (zero weight)
        undefined, // missing kappa (zero weight)
        f1_only_config,
      )

      // Should be equal to the F1 value
      expect(f1_only_score).toBe(0.8)

      // Test with only RMSD available in RMSD-only config
      const rmsd_only_score = calculate_combined_score(
        undefined, // missing F1 (zero weight)
        0.05, // good RMSD
        undefined, // missing kappa (zero weight)
        rmsd_only_config,
      )

      // RMSD is inverted and normalized to [0,1]
      expect(rmsd_only_score).toBeGreaterThan(0.7)

      // Test with only kappa available in kappa-only config
      const kappa_only_score = calculate_combined_score(
        undefined, // missing F1 (zero weight)
        undefined, // missing RMSD (zero weight)
        0.3, // good kappa
        kappa_only_config,
      )

      // Kappa is normalized to [0,1]
      expect(kappa_only_score).toBeGreaterThan(0.7)
    })

    it(`returns NaN when weighted metrics are missing`, () => {
      // Create a config with non-zero weights
      const config = {
        ...DEFAULT_COMBINED_METRIC_CONFIG,
        weights: [
          { metric: `F1`, value: 1, display: `F1`, description: `` },
          { metric: `RMSD`, value: 0, display: `RMSD`, description: `` },
          { metric: `kappa_SRME`, value: 0, display: `kappa`, description: `` },
        ],
      }

      // Missing F1 but F1 weight is 1
      const score = calculate_combined_score(undefined, 0.05, 0.3, config)

      expect(isNaN(score)).toBe(true)
    })

    it(`correctly weights metrics according to config`, () => {
      const f1 = 1.0 // perfect F1
      const rmsd = 0.3 // poor RMSD (maximum baseline value)
      const kappa = 2.0 // poor kappa (maximum value)

      // Test with equal weights
      const equal_weights = {
        ...DEFAULT_COMBINED_METRIC_CONFIG,
        weights: [
          { metric: `F1`, value: 1 / 3, display: `F1`, description: `` },
          { metric: `RMSD`, value: 1 / 3, display: `RMSD`, description: `` },
          { metric: `kappa_SRME`, value: 1 / 3, display: `kappa`, description: `` },
        ],
      }

      const equal_score = calculate_combined_score(f1, rmsd, kappa, equal_weights)

      // Perfect F1 (1.0), worst RMSD (0.0), worst kappa (0.0)
      // Equal weights: (1.0 + 0.0 + 0.0) / 3 = 0.333...
      expect(equal_score).toBeCloseTo(1 / 3, 1)

      // Test with F1-only weight
      const f1_only_weights = {
        ...DEFAULT_COMBINED_METRIC_CONFIG,
        weights: [
          { metric: `F1`, value: 1, display: `F1`, description: `` },
          { metric: `RMSD`, value: 0, display: `RMSD`, description: `` },
          { metric: `kappa_SRME`, value: 0, display: `kappa`, description: `` },
        ],
      }

      const f1_weighted_score = calculate_combined_score(f1, rmsd, kappa, f1_only_weights)

      // Should be equal to F1 value (1.0)
      expect(f1_weighted_score).toBe(1.0)
    })

    it(`normalizes RMSD correctly`, () => {
      // Create a config with only RMSD weighted
      const rmsd_only_config = {
        ...DEFAULT_COMBINED_METRIC_CONFIG,
        weights: [
          { metric: `F1`, value: 0, display: `F1`, description: `` },
          { metric: `RMSD`, value: 1, display: `RMSD`, description: `` },
          { metric: `kappa_SRME`, value: 0, display: `kappa`, description: `` },
        ],
      }

      // Test with excellent RMSD (close to 0)
      const excellent_score = calculate_combined_score(
        undefined,
        0.01,
        undefined,
        rmsd_only_config,
      )
      expect(excellent_score).toBeGreaterThan(0.9)

      // Test with poor RMSD (at baseline)
      const poor_score = calculate_combined_score(
        undefined,
        0.3,
        undefined,
        rmsd_only_config,
      )
      expect(poor_score).toBeCloseTo(0, 1)

      // Test with mid-range RMSD
      const mid_score = calculate_combined_score(
        undefined,
        0.15,
        undefined,
        rmsd_only_config,
      )
      expect(mid_score).toBeCloseTo(0.5, 1)
    })

    it(`normalizes kappa SRME correctly`, () => {
      // Create a config with only kappa weighted
      const kappa_only_config = {
        ...DEFAULT_COMBINED_METRIC_CONFIG,
        weights: [
          { metric: `F1`, value: 0, display: `F1`, description: `` },
          { metric: `RMSD`, value: 0, display: `RMSD`, description: `` },
          { metric: `kappa_SRME`, value: 1, display: `kappa`, description: `` },
        ],
      }

      // Test with excellent kappa (close to 0)
      const excellent_score = calculate_combined_score(
        undefined,
        undefined,
        0.1,
        kappa_only_config,
      )
      expect(excellent_score).toBeGreaterThan(0.9)

      // Test with poor kappa (maximum value is 2)
      const poor_score = calculate_combined_score(
        undefined,
        undefined,
        2.0,
        kappa_only_config,
      )
      expect(poor_score).toBe(0)

      // Test with mid-range kappa
      const mid_score = calculate_combined_score(
        undefined,
        undefined,
        1.0,
        kappa_only_config,
      )
      expect(mid_score).toBeCloseTo(0.5, 1)
    })

    it(`assigns correct default weights`, () => {
      expect(F1_DEFAULT_WEIGHT).toBe(0.5)
      expect(RMSD_DEFAULT_WEIGHT).toBe(0.1)
      expect(KAPPA_DEFAULT_WEIGHT).toBe(0.4)

      const sum = F1_DEFAULT_WEIGHT + RMSD_DEFAULT_WEIGHT + KAPPA_DEFAULT_WEIGHT
      expect(sum).toBe(1.0)
    })
  })
})
