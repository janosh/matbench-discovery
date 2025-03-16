import { RadarChart } from '$lib'
import type { MetricWeight } from '$lib/types'
import { mount, tick } from 'svelte'
import { describe, expect, it, vi } from 'vitest'

describe(`RadarChart`, () => {
  // Model metrics data
  const model_metrics = {
    test_model: {
      f1: 0.85,
      rmsd: 0.12,
      kappa_SRME: 0.28,
    },
  }

  // Define weights as MetricWeight[] as expected by the component
  const sample_weights: MetricWeight[] = [
    {
      metric: `f1`,
      label: `F1`,
      description: `F1 score for model performance`,
      value: 0.5,
    },
    {
      metric: `rmsd`,
      label: `RMSD`,
      description: `Root mean square displacement`,
      value: 0.3,
    },
    {
      metric: `kappa_SRME`,
      label: `κ<sub>SRME</sub>`,
      description: `Symmetric relative mean error for thermal conductivity`,
      value: 0.2,
    },
  ]

  it.skip(`handles weight changes`, async () => {
    // Note: Skipping this test as JSDOM doesn't properly handle SVG DOM events
    // This would need to be tested in a browser environment with proper event simulation

    const onChange = vi.fn()

    mount(RadarChart, {
      target: document.body,
      props: {
        metrics: model_metrics,
        weights: sample_weights,
        onChange,
      },
    })

    await tick()

    // In a real browser environment, we would:
    // 1. Find a weight knob
    // 2. Simulate drag events
    // 3. Verify onChange is called with updated weights

    // For now, we just check that the knobs are rendered
    const knobs = document.body.querySelectorAll(`.weight-knob`)
    expect(knobs.length).toBe(3)
  })

  it.skip(`handles dragging interaction`, async () => {
    // Note: Skipping this test as JSDOM has limitations with SVG+mouse events
    // This would need to be tested in a browser environment

    const onChange = vi.fn()

    mount(RadarChart, {
      target: document.body,
      props: {
        metrics: model_metrics,
        weights: sample_weights,
        onChange,
      },
    })

    await tick()

    // In a real browser environment, we would:
    // 1. Find a weight knob
    // 2. Simulate mousedown, mousemove, mouseup events
    // 3. Verify the knob position changes and onChange is called

    // For now, we just check that the chart is interactive
    const svg = document.body.querySelector(`svg`)
    expect(svg).toBeDefined()
  })

  it(`renders with custom size`, async () => {
    const size = 400

    mount(RadarChart, {
      target: document.body,
      props: {
        metrics: model_metrics,
        weights: sample_weights,
        size,
      },
    })

    await tick()

    const svg = document.body.querySelector(`svg`)
    expect(svg).toBeDefined()
    expect(svg?.getAttribute(`width`)).toBe(String(size))
    expect(svg?.getAttribute(`height`)).toBe(String(size))
  })
})
