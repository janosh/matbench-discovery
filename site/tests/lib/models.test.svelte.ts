import { CPS_CONFIG } from '$lib/combined_perf_score.svelte'
import {
  calculate_training_sizes,
  MODEL_METADATA_PATHS,
  MODELS,
  update_models_cps,
} from '$lib/models.svelte'
import { describe, expect, it, vi } from 'vitest'

describe(`calculate_training_sizes`, () => {
  it(`should return zeros for empty training sets`, () => {
    const result = calculate_training_sizes([])
    expect(result.total_materials).toBe(0)
    expect(result.total_structures).toBe(0)
  })

  it(`should calculate totals correctly for a single dataset`, () => {
    const result = calculate_training_sizes([`MP 2022`])
    expect(result.total_materials).toBe(154_719)
    expect(result.total_structures).toBe(154_719)
  })

  it(`should calculate totals correctly for multiple datasets`, () => {
    const result = calculate_training_sizes([`MP 2022`, `MPtrj`])
    expect(result.total_materials).toBe(300_642)
    expect(result.total_structures).toBe(1_735_114)
  })

  it(`should use n_structures as n_materials when n_materials is not specified`, () => {
    const result = calculate_training_sizes([`MPF`])
    expect(result.total_materials).toBe(62_783) // Should equal n_structures
    expect(result.total_structures).toBe(188_349)
  })

  it(`should skip datasets that don't exist and log a warning`, () => {
    const console_warn_spy = vi.spyOn(console, `warn`).mockImplementation(() => {})

    const result = calculate_training_sizes([`NonExistentDataset`, `MP 2022`])

    expect(console_warn_spy).toHaveBeenCalledWith(
      `Training set NonExistentDataset not found in DATASETS`,
    )
    expect(result.total_materials).toBe(154_719) // Only from MP 2022
    expect(result.total_structures).toBe(154_719) // Only from MP 2022

    console_warn_spy.mockRestore()
  })
})

describe(`MODELS array`, () => {
  it(`should be defined and be an array`, () => {
    expect(MODELS).toBeDefined()
    expect(Array.isArray(MODELS)).toBe(true)
  })

  it(`should have processed models with calculated properties`, () => {
    // Skip if no models available in test environment
    if (MODELS.length === 0) return

    // Check that models have the expected structure
    const model = MODELS[0]
    expect(model).toHaveProperty(`dirname`)
    expect(model).toHaveProperty(`metadata_file`)
    expect(model).toHaveProperty(`color`)
    expect(model).toHaveProperty(`n_training_materials`)
    expect(model).toHaveProperty(`n_training_structures`)
  })
})

describe(`MODEL_METADATA_PATHS`, () => {
  it(`should be defined and be an object`, () => {
    expect(MODEL_METADATA_PATHS).toBeDefined()
    expect(typeof MODEL_METADATA_PATHS).toBe(`object`)
  })
})

describe(`CPS_CONFIG`, () => {
  it(`should be defined and match DEFAULT_CPS_CONFIG initially`, () => {
    expect(CPS_CONFIG).toBeDefined()
    expect(CPS_CONFIG.F1.weight).toBeDefined()
    expect(CPS_CONFIG.RMSD.weight).toBeDefined()
    expect(CPS_CONFIG.κ_SRME.weight).toBeDefined()
  })

  it(`should be reactive (modifiable)`, () => {
    // Store original weights
    const original_f1_weight = CPS_CONFIG.F1.weight

    // Modify weight
    CPS_CONFIG.F1.weight = 0.8

    // Check that it was updated
    expect(CPS_CONFIG.F1.weight).toBe(0.8)

    // Restore original weight
    CPS_CONFIG.F1.weight = original_f1_weight
  })
})

describe(`update_models_cps`, () => {
  it(`should update CPS for models based on metrics and current weights`, () => {
    // Skip test if no models available
    if (MODELS.length === 0) return

    // Act: Call the function under test
    update_models_cps(MODELS, CPS_CONFIG)

    // Assert: Check if at least one model has a non-NaN CPS value
    const models_with_cps = MODELS.filter((model) => !isNaN(Number(model.CPS)))
    expect(models_with_cps.length).toBeGreaterThan(0)
  })

  it(`should set CPS to NaN when required metrics are missing`, () => {
    // Skip test if no models available
    if (MODELS.length === 0) return

    // Act: Set weights that would require all metrics to be present
    CPS_CONFIG.F1.weight = 0.3
    CPS_CONFIG.RMSD.weight = 0.3
    CPS_CONFIG.κ_SRME.weight = 0.4
    update_models_cps(MODELS, CPS_CONFIG)

    // Assert: Check if some models have NaN CPS (due to missing metrics)
    const models_with_nan_cps = MODELS.filter((model) => isNaN(Number(model.CPS)))
    expect(models_with_nan_cps.length).toBeGreaterThan(0)
  })

  it(`should handle different weight configurations correctly`, () => {
    // First configuration - only F1 matters
    CPS_CONFIG.F1.weight = 1.0
    CPS_CONFIG.RMSD.weight = 0.0
    CPS_CONFIG.κ_SRME.weight = 0.0
    update_models_cps(MODELS, CPS_CONFIG)

    const f1_cps_values = MODELS.map((model) => Number(model.CPS))

    // Second configuration - only RMSD matters
    CPS_CONFIG.F1.weight = 0.0
    CPS_CONFIG.RMSD.weight = 1.0
    CPS_CONFIG.κ_SRME.weight = 0.0
    update_models_cps(MODELS, CPS_CONFIG)

    const rmsd_cps_values = MODELS.map((model) => Number(model.CPS))

    expect(f1_cps_values).not.toEqual(rmsd_cps_values)

    CPS_CONFIG.F1.weight = 0.5
    CPS_CONFIG.RMSD.weight = 0
    CPS_CONFIG.κ_SRME.weight = 0.5
    update_models_cps(MODELS, CPS_CONFIG)

    const f1_rmsd_cps_values = MODELS.map((model) => Number(model.CPS))
    expect(f1_rmsd_cps_values).not.toEqual(f1_cps_values)
    expect(f1_rmsd_cps_values).not.toEqual(rmsd_cps_values)
  })
})
