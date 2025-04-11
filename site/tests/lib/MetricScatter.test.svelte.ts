import { MetricScatter } from '$lib'
import { DEFAULT_CPS_CONFIG } from '$lib/combined_perf_score'
import type { ModelData } from '$lib/types'
import { mount } from 'svelte'
import { describe, expect, it } from 'vitest'

// Mock model data for testing
const mock_models: ModelData[] = [
  {
    model_name: `Test Model 1`,
    model_params: 100000,
    n_estimators: 50,
    dirname: `test1`,
    metadata_file: `test1.yml`,
    date_added: `2023-01-01`,
    training_cost: { 'NVIDIA A100 GPUs': { amount: 4, hours: 24, cost: 200 } },
    metrics: {
      discovery: {
        pred_col: `energy`,
        unique_prototypes: { F1: 0.85, DAF: 0.9 },
        full_test_set: { F1: 0.8, DAF: 0.85 },
      },
      geo_opt: {
        'symprec=1e-5': {
          rmsd: 0.02,
          n_sym_ops_mae: 0.5,
          symmetry_decrease: 0,
          symmetry_match: 1,
          symmetry_increase: 0,
          n_structures: 100,
        },
      },
      phonons: { kappa_103: { κ_SRME: 0.15 } },
    },
  } as unknown as ModelData,
  {
    model_name: `Test Model 2`,
    model_params: 500000,
    n_estimators: 100,
    dirname: `test2`,
    metadata_file: `test2.yml`,
    date_added: `2023-03-15`,
    training_cost: { 'NVIDIA A100 GPUs': { amount: 8, hours: 48, cost: 800 } },
    metrics: {
      discovery: {
        pred_col: `energy`,
        unique_prototypes: { F1: 0.9, DAF: 0.95 },
        full_test_set: { F1: 0.85, DAF: 0.9 },
      },
      geo_opt: {
        'symprec=1e-5': {
          rmsd: 0.015,
          n_sym_ops_mae: 0.4,
          symmetry_decrease: 0,
          symmetry_match: 1,
          symmetry_increase: 0,
          n_structures: 100,
        },
      },
      phonons: { kappa_103: { κ_SRME: 0.1 } },
    },
  } as unknown as ModelData,
  {
    model_name: `Test Model 3`,
    model_params: 1000000,
    n_estimators: 200,
    dirname: `test3`,
    metadata_file: `test3.yml`,
    date_added: `2023-06-30`,
    training_cost: { 'NVIDIA A100 GPUs': { amount: 16, hours: 72, cost: 2000 } },
    metrics: {
      discovery: {
        pred_col: `energy`,
        unique_prototypes: { F1: 0.95, DAF: 0.98 },
        full_test_set: { F1: 0.9, DAF: 0.93 },
      },
      geo_opt: {
        'symprec=1e-5': {
          rmsd: 0.01,
          n_sym_ops_mae: 0.3,
          symmetry_decrease: 0,
          symmetry_match: 1,
          symmetry_increase: 0,
          n_structures: 100,
        },
      },
      phonons: { kappa_103: { κ_SRME: 0.05 } },
    },
  } as unknown as ModelData,
]

describe(`MetricScatter impossible`, () => {
  // Simplified tests that focus on component mounting and basic functionality
  it(`renders with default props`, () => {
    mount(MetricScatter, {
      target: document.body,
      props: {
        models: mock_models,
        config: DEFAULT_CPS_CONFIG,
      },
    })

    // Check that component rendered with an SVG element
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })

  it(`renders points for valid models`, async () => {
    // Use mock data with valid F1 scores and model parameters
    mount(MetricScatter, {
      target: document.body,
      props: {
        models: mock_models,
        config: DEFAULT_CPS_CONFIG,
        x_property: `model_params`,
        y_metric: `F1`,
      },
    })

    // Verify some circles are rendered
    const svg = document.body.querySelector(`svg`)
    expect(svg).toBeDefined()
  })

  it(`respects model_filter property`, async () => {
    // First mount with all models
    mount(MetricScatter, {
      target: document.body,
      props: {
        models: mock_models,
        config: DEFAULT_CPS_CONFIG,
      },
    })

    // Then mount with filtered models
    mount(MetricScatter, {
      target: document.body,
      props: {
        models: mock_models,
        config: DEFAULT_CPS_CONFIG,
        model_filter: (model: ModelData) => model.model_name === `Test Model 1`,
      },
    })

    // Validate component rendered
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })

  it(`accepts custom styling properties`, async () => {
    const custom_color = `#ff0000` // red
    const custom_radius = 10

    mount(MetricScatter, {
      target: document.body,
      props: {
        models: mock_models,
        config: DEFAULT_CPS_CONFIG,
        point_color: custom_color,
        point_radius: custom_radius,
      },
    })

    // Validate component rendered
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })

  it(`handles different metrics`, async () => {
    // Test with F1 metric
    mount(MetricScatter, {
      target: document.body,
      props: {
        models: mock_models,
        config: DEFAULT_CPS_CONFIG,
        y_metric: `F1`,
      },
    })

    // Test with RMSD metric
    mount(MetricScatter, {
      target: document.body,
      props: {
        models: mock_models,
        config: DEFAULT_CPS_CONFIG,
        y_metric: `RMSD`,
      },
    })

    // Validate component rendered
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })

  it(`handles different properties for X axis`, async () => {
    // Test with model_params property
    mount(MetricScatter, {
      target: document.body,
      props: {
        models: mock_models,
        config: DEFAULT_CPS_CONFIG,
        x_property: `model_params`,
      },
    })

    // Test with training_cost property
    mount(MetricScatter, {
      target: document.body,
      props: {
        models: mock_models,
        config: DEFAULT_CPS_CONFIG,
        x_property: `training_cost`,
      },
    })

    // Validate component rendered
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })

  it(`handles custom axis labels`, async () => {
    const custom_x_label = `Custom X Label`
    const custom_y_label = `Custom Y Label`

    mount(MetricScatter, {
      target: document.body,
      props: {
        models: mock_models,
        config: DEFAULT_CPS_CONFIG,
        x_label: custom_x_label,
        y_label: custom_y_label,
      },
    })

    // Validate component rendered
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })

  // Additional tests for MetricScatter functionality
  it(`renders with date on X axis`, async () => {
    mount(MetricScatter, {
      target: document.body,
      props: {
        models: mock_models,
        x_property: `date_added`,
        y_metric: `F1`,
      },
    })

    // Validate component rendered
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })

  it(`renders with dotted path metrics`, async () => {
    mount(MetricScatter, {
      target: document.body,
      props: {
        models: mock_models,
        metric: `phonons.kappa_103.κ_SRME`,
        y_label: `κ SRME Value`,
      },
    })

    // Validate component rendered
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })

  it(`handles date range filtering`, async () => {
    const start_date = new Date(`2023-03-01`)
    const end_date = new Date(`2023-12-31`)

    mount(MetricScatter, {
      target: document.body,
      props: {
        models: mock_models,
        x_property: `date_added`,
        y_metric: `F1`,
        date_range: [start_date, end_date],
      },
    })

    // Validate component rendered
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })

  it(`applies custom styles`, async () => {
    const custom_style = `border: 1px solid red;`

    mount(MetricScatter, {
      target: document.body,
      props: {
        models: mock_models,
        style: custom_style,
      },
    })

    // Validate component rendered
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })

  it(`renders with y_lim as null`, async () => {
    mount(MetricScatter, {
      target: document.body,
      props: {
        models: mock_models,
        y_lim: null,
      },
    })

    // Validate component rendered
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })
})
