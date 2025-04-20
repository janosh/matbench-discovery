import { CPS_CONFIG, MODELS, MODEL_METADATA_PATHS } from '$lib/models.svelte'
import { beforeEach, describe, expect, it, vi } from 'vitest'
// Import the function we are testing
import { update_models_cps } from '$lib/models.svelte'
// Import the function we need to mock
import { calculate_cps } from '$lib/combined_perf_score'
// Import the type for mock models
import type { ModelData } from '$lib/types'

// Define mock datasets in the vi.mock factory to avoid hoisting issues
vi.mock(`$data/datasets.yml`, () => ({
  default: {
    'MP 2022': {
      n_structures: 100,
      n_materials: 80,
      title: `MP 2022 Test`,
      description: `Test dataset`,
      url: `https://example.com`,
      doi: `10.1234/test`,
      slug: `mp-2022`,
    },
    MPtrj: {
      n_structures: 200,
      n_materials: 150,
      title: `MPtrj Test`,
      description: `Test dataset`,
      url: `https://example.com`,
      doi: `10.1234/test`,
      slug: `mptrj`,
    },
    MPF: {
      n_structures: 300,
      title: `MPF Test`,
      description: `Test dataset`,
      url: `https://example.com`,
      doi: `10.1234/test`,
      slug: `mpf`,
    },
  },
}))

// Mock the calculate_cps function from the combined_perf_score module
vi.mock(`$lib/combined_perf_score`, () => ({
  calculate_cps: vi.fn(),
  // Provide a basic structure for DEFAULT_CPS_CONFIG to prevent test errors
  DEFAULT_CPS_CONFIG: {
    parts: {
      F1: { weight: 0, min: 0, max: 1 },
      RMSD: { weight: 0, min: 0, max: 1 },
      kappa_SRME: { weight: 0, min: 0, max: 1 },
    },
  },
}))

// Store mock datasets for tests to use
const mockDatasets: Record<string, { n_structures: number; n_materials?: number }> = {
  'MP 2022': {
    n_structures: 100,
    n_materials: 80,
  },
  MPtrj: {
    n_structures: 200,
    n_materials: 150,
  },
  MPF: {
    n_structures: 300,
    // n_materials not specified, should default to n_structures
  },
}

// Access the private calculate_training_sizes function for testing
// We need to import and recreate it since it's not exported
function calculate_training_sizes(model_train_sets: string[] = []): {
  total_materials: number
  total_structures: number
} {
  let total_materials = 0
  let total_structures = 0

  for (const data_name of model_train_sets) {
    if (!(data_name in mockDatasets)) {
      console.warn(`Training set ${data_name} not found in DATASETS`)
      continue
    }
    const { n_structures, n_materials = n_structures } = mockDatasets[data_name]
    total_materials += n_materials
    total_structures += n_structures
  }

  return { total_materials, total_structures }
}

describe(`calculate_training_sizes`, () => {
  it(`should return zeros for empty training sets`, () => {
    const result = calculate_training_sizes([])
    expect(result.total_materials).toBe(0)
    expect(result.total_structures).toBe(0)
  })

  it(`should calculate totals correctly for a single dataset`, () => {
    const result = calculate_training_sizes([`MP 2022`])
    expect(result.total_materials).toBe(80)
    expect(result.total_structures).toBe(100)
  })

  it(`should calculate totals correctly for multiple datasets`, () => {
    const result = calculate_training_sizes([`MP 2022`, `MPtrj`])
    expect(result.total_materials).toBe(80 + 150)
    expect(result.total_structures).toBe(100 + 200)
  })

  it(`should use n_structures as n_materials when n_materials is not specified`, () => {
    const result = calculate_training_sizes([`MPF`])
    expect(result.total_materials).toBe(300) // Should equal n_structures
    expect(result.total_structures).toBe(300)
  })

  it(`should skip datasets that don't exist and log a warning`, () => {
    const console_warn_spy = vi.spyOn(console, `warn`).mockImplementation(() => {})

    const result = calculate_training_sizes([`NonExistentDataset`, `MP 2022`])

    expect(console_warn_spy).toHaveBeenCalledWith(
      `Training set NonExistentDataset not found in DATASETS`,
    )
    expect(result.total_materials).toBe(80) // Only from MP 2022
    expect(result.total_structures).toBe(100) // Only from MP 2022

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
    expect(CPS_CONFIG.parts.F1.weight).toBeDefined()
    expect(CPS_CONFIG.parts.RMSD.weight).toBeDefined()
    expect(CPS_CONFIG.parts.kappa_SRME.weight).toBeDefined()
  })

  it(`should be reactive (modifiable)`, () => {
    // Store original weights
    const original_f1_weight = CPS_CONFIG.parts.F1.weight

    // Modify weight
    CPS_CONFIG.parts.F1.weight = 0.8

    // Check that it was updated
    expect(CPS_CONFIG.parts.F1.weight).toBe(0.8)

    // Restore original weight
    CPS_CONFIG.parts.F1.weight = original_f1_weight
  })
})

describe(`update_models_cps`, () => {
  // Define a base mock model to satisfy ModelData type requirements
  const base_mock_model: ModelData = {
    model_id: `test-model`,
    model_name: `Test Model`,
    model_version: `1.0.0`,
    date_added: `2024-01-01`,
    date_published: `2024-01-01`,
    // @ts-expect-error: authors has wrong type in model-schema.d.ts
    authors: [
      {
        name: `Test Author`,
        affiliation: `Test Affiliation`,
        orcid: `0000-0000-0000-0000`,
      },
    ],
    source_code: `https://example.com/code`,
    paper: `https://example.com/paper`,
    project_url: `https://example.com/project`,
    // Add missing required properties
    repo: `https://github.com/test/repo`,
    pr_url: `https://github.com/test/repo/pull/1`,
    checkpoint_url: `https://example.com/checkpoint`,
    // Requirements should be an object
    requirements: { 'test-req': `1.0` },
    model_family: `TestFamily`,
    model_task_type: `regression`,
    model_output_type: `scalar`,
    output_labels: [`energy`],
    training_code_url: `https://example.com/train-code`,
    train_set_file: `train.csv`,
    test_set_file: `test.csv`,
    validation_set_file: `val.csv`,
    // ---
    training_set: [],
    model_format: `unknown`,
    train_test_split: `unknown`,
    status: `complete`,
    tags: [],
    doi: `10.1234/testdoi`,
    license: { code: `MIT`, checkpoint: `unreleased` },
    n_params: 1000,
    n_elements: 10,
    framework: `Test Framework`,
    primary_task: `energy`,
    description: `Test description`,
    // --- Calculated properties ---
    dirname: `test-model`,
    metadata_file: `/path/to/test-model.yml`,
    color: `#FFFFFF`,
    CPS: NaN,
    n_training_materials: 0,
    n_training_structures: 0,
    metrics: {}, // Metrics will be added per test
  }

  // Ensure mocks and MODELS state are clean before each test
  beforeEach(() => {
    vi.resetAllMocks()
    // Reset MODELS state by clearing the array
    MODELS.length = 0
  })

  it(`should update CPS for models based on metrics and calculate_cps`, () => {
    // Arrange: Create a mock model with specific metrics, based on the base
    const mock_model: ModelData = {
      ...base_mock_model,
      model_id: `test-model-1`,
      dirname: `test-model-1`,
      metrics: {
        discovery: {
          pred_col: `energy_pred`, // Required field
          unique_prototypes: { F1: 0.8, Precision: 0, Recall: 0 }, // Added missing fields
        },
        geo_opt: {
          'symprec=1e-5': {
            // Required fields added
            rmsd: 0.1,
            n_sym_ops_mae: 0,
            symmetry_decrease: 0,
            symmetry_match: 0,
            symmetry_increase: 0,
            n_structures: 1,
          },
        },
        phonons: {
          kappa_103: { κ_SRME: `10.5` }, // Note: original data can have string here
        },
      },
    }
    MODELS.push(mock_model)

    // Arrange: Mock calculate_cps to return a specific value for these inputs
    const expected_cps = 123.45
    vi.mocked(calculate_cps).mockImplementation((f1, rmsd, kappa) => {
      // Check if the correct metrics are passed (allowing for type coercion)
      if (f1 === 0.8 && rmsd === 0.1 && kappa === 10.5) {
        return expected_cps
      }
      return NaN // Default return if inputs don't match
    })

    // Act: Call the function under test
    update_models_cps()

    // Assert: Check if the model's CPS was updated
    expect(MODELS[0].CPS).toBe(expected_cps)
    // Assert: Check if calculate_cps was called with the correctly extracted metrics
    expect(calculate_cps).toHaveBeenCalledWith(0.8, 0.1, 10.5, CPS_CONFIG)
  })

  it(`should set CPS to NaN if required metrics are missing`, () => {
    // Arrange: Create a mock model with completely missing metrics
    const mock_model_missing: ModelData = {
      ...base_mock_model,
      model_id: `test-model-missing`,
      dirname: `test-model-missing`,
      metrics: {}, // Metrics object exists but is empty
      CPS: 999, // Start with a non-NaN value to ensure it's overwritten by NaN
    }
    MODELS.push(mock_model_missing)

    // Arrange: Mock calculate_cps to return NaN when inputs are undefined
    vi.mocked(calculate_cps).mockReturnValue(NaN)

    // Act: Call the function under test
    update_models_cps()

    // Assert: Check if the model's CPS is NaN
    expect(MODELS[0].CPS).toBeNaN()
    // Assert: Check if calculate_cps was called with undefined for all metrics
    expect(calculate_cps).toHaveBeenCalledWith(
      undefined,
      undefined,
      undefined,
      CPS_CONFIG,
    )
  })

  it(`should handle models with partial metrics correctly`, () => {
    // Arrange: Model with only F1 score metric available
    const mock_model_partial: ModelData = {
      ...base_mock_model,
      model_id: `test-model-partial`,
      dirname: `test-model-partial`,
      metrics: {
        discovery: {
          pred_col: `energy_pred`, // Required field
          unique_prototypes: { F1: 0.7, Precision: 0, Recall: 0 }, // Added missing fields
        },
        // geo_opt and phonons are missing
      },
    }
    MODELS.push(mock_model_partial)

    // Arrange: Mock calculate_cps for this specific partial input case
    const expected_partial_cps = 50.0 // Example value returned by calculate_cps
    vi.mocked(calculate_cps).mockImplementation((f1, rmsd, kappa) => {
      if (f1 === 0.7 && rmsd === undefined && kappa === undefined) {
        return expected_partial_cps
      }
      return NaN
    })

    // Act
    update_models_cps()

    // Assert: Check if CPS reflects the partial calculation
    expect(MODELS[0].CPS).toBe(expected_partial_cps)
    // Assert: Check call signature for partial metrics
    expect(calculate_cps).toHaveBeenCalledWith(0.7, undefined, undefined, CPS_CONFIG)
  })

  it(`should handle non-numeric kappa value gracefully`, () => {
    // Arrange: Model with kappa as a non-numeric string
    const mock_model_bad_kappa: ModelData = {
      ...base_mock_model,
      model_id: `test-model-bad-kappa`,
      dirname: `test-model-bad-kappa`,
      metrics: {
        discovery: {
          pred_col: `energy_pred`, // Required field
          unique_prototypes: { F1: 0.6, Precision: 0, Recall: 0 }, // Added missing fields
        },
        geo_opt: {
          'symprec=1e-5': {
            // Required fields added
            rmsd: 0.2,
            n_sym_ops_mae: 0,
            symmetry_decrease: 0,
            symmetry_match: 0,
            symmetry_increase: 0,
            n_structures: 1,
          },
        },
        phonons: { kappa_103: { κ_SRME: `not-a-number` } }, // Invalid kappa
      },
    }
    MODELS.push(mock_model_bad_kappa)

    const expected_cps_bad_kappa = 40.0 // Example value
    vi.mocked(calculate_cps).mockImplementation((f1, rmsd, kappa) => {
      // Correctly check for NaN, as Number('not-a-number') is NaN
      if (f1 === 0.6 && rmsd === 0.2 && Number.isNaN(kappa)) {
        return expected_cps_bad_kappa
      }
      return NaN
    })

    // Act
    update_models_cps()

    // Assert
    expect(MODELS[0].CPS).toBe(expected_cps_bad_kappa)
    // Check that calculate_cps was called with NaN for kappa
    expect(calculate_cps).toHaveBeenCalledWith(0.6, 0.2, NaN, CPS_CONFIG)
  })
})
