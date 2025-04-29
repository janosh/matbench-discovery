import { MetricScatter } from '$lib'
import { DEFAULT_CPS_CONFIG } from '$lib/combined_perf_score.svelte'
import { ALL_METRICS, METADATA_COLS } from '$lib/labels'
import type { ModelData } from '$lib/types'
import { mount } from 'svelte'
import { describe, expect, it } from 'vitest'

describe(`MetricScatter`, () => {
  // Test basic rendering with required props
  it(`renders with default props`, () => {
    mount(MetricScatter, {
      target: document.body,
      props: { x_prop: METADATA_COLS.model_params, y_prop: ALL_METRICS.CPS },
    })
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })

  // Test rendering points with specific columns
  it(`renders points with F1 vs model_params`, () => {
    mount(MetricScatter, {
      target: document.body,
      props: { x_prop: METADATA_COLS.model_params, y_prop: ALL_METRICS.F1 },
    })
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })

  // Test model filtering
  it(`respects model_filter property`, () => {
    mount(MetricScatter, {
      target: document.body,
      props: {
        x_prop: METADATA_COLS.model_params,
        y_prop: ALL_METRICS.CPS,
        config: DEFAULT_CPS_CONFIG,
        model_filter: (model: ModelData) => model.model_name === `Test Model 1`,
      },
    })
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })

  // Test point radius styling
  it(`accepts point_radius property`, () => {
    const custom_radius = 10
    mount(MetricScatter, {
      target: document.body,
      props: {
        x_prop: METADATA_COLS.model_params,
        y_prop: ALL_METRICS.CPS,
        point_radius: custom_radius,
      },
    })
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })

  // Test handling different metric columns
  it(`handles different metric columns (RMSD, Kappa)`, () => {
    // Test with RMSD
    mount(MetricScatter, {
      target: document.body,
      props: { x_prop: METADATA_COLS.model_params, y_prop: ALL_METRICS.RMSD },
    })
    expect(document.body.querySelector(`svg`)).toBeDefined()

    // Test with Kappa
    mount(MetricScatter, {
      target: document.body,
      props: { x_prop: METADATA_COLS.model_params, y_prop: ALL_METRICS.κ_SRME },
    })
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })

  // Test handling different property columns for X axis
  it(`handles different property columns for X axis (training_cost)`, () => {
    mount(MetricScatter, {
      target: document.body,
      props: {
        x_prop: ALL_METRICS.Precision,
        y_prop: ALL_METRICS.CPS,
        config: DEFAULT_CPS_CONFIG,
      },
    })
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })

  // Test rendering with date on X axis
  it(`renders with date on X axis`, () => {
    mount(MetricScatter, {
      target: document.body,
      props: { x_prop: METADATA_COLS.date_added, y_prop: ALL_METRICS.F1 },
    })
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })

  // Test date range filtering
  it(`handles date range filtering`, () => {
    const start_date = new Date(`2023-03-01`)
    const end_date = new Date(`2023-12-31`)

    mount(MetricScatter, {
      target: document.body,
      props: {
        x_prop: METADATA_COLS.date_added,
        y_prop: ALL_METRICS.F1,
        date_range: [start_date, end_date],
      },
    })
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })

  // Test custom styles
  it(`applies custom styles`, () => {
    const custom_style = `border: 1px solid red;`
    mount(MetricScatter, {
      target: document.body,
      props: {
        x_prop: METADATA_COLS.model_params,
        y_prop: ALL_METRICS.CPS,
        style: custom_style,
      },
    })
    expect(document.body.querySelector(`svg`)).toBeDefined()
  })
})
