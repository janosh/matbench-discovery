import { RadarChart } from '$lib'
import { CPS_CONFIG, DEFAULT_CPS_CONFIG } from '$lib/combined_perf_score.svelte'
import { ALL_METRICS } from '$lib/labels'
import { update_models_cps } from '$lib/models.svelte'
import { mount } from 'svelte'
import { beforeEach, describe, expect, it, vi } from 'vitest'

describe(`RadarChart`, () => {
  beforeEach(() => {
    // Mock the update_models_cps function to avoid side effects
    vi.mock(`$lib/models.svelte`, async () => {
      const actual = await vi.importActual(`$lib/models.svelte`)
      return { ...actual, update_models_cps: vi.fn() }
    })
  })

  it(`renders with default props`, () => {
    mount(RadarChart, { target: document.body })

    // Check that the component rendered with an SVG element
    const svg = document.body.querySelector(`svg`)
    expect(svg).toBeDefined()

    // Check that all three axis lines are rendered (for F1, kappa, RMSD)
    const axis_lines = document.body.querySelectorAll(`line`)
    expect(axis_lines).toHaveLength(3)
    // Verify stroke properties for each axis line
    axis_lines.forEach((line) => {
      expect(line.getAttribute(`stroke`)).toBe(`var(--border)`)
      expect(line.getAttribute(`stroke-width`)).toBe(`1`)
    })

    // Check that the triangle area is rendered
    const triangle_area = document.body.querySelector(
      `path[fill="var(--nav-bg)"]`,
    )
    expect(triangle_area).toBeDefined()

    // Check that the draggable point is rendered
    const draggable_point = document.body.querySelector(`circle[role="button"]`)
    expect(draggable_point).toBeDefined()
  })

  it(`accepts size prop`, () => {
    const custom_size = 300

    mount(RadarChart, { target: document.body, props: { size: custom_size } })

    const svg = document.body.querySelector(`svg`)
    expect(svg).toBeDefined()

    // Verify that at least one circle is rendered, which indirectly confirms
    // that the component initialized properly with the size prop
    const circles = document.body.querySelectorAll(`circle`)
    expect(circles.length).toBeGreaterThan(0)
  })

  it(`resets weights when reset button is clicked`, () => {
    mount(RadarChart, { target: document.body })

    // Modify weights to non-default values
    CPS_CONFIG.F1.weight = 0.6
    CPS_CONFIG.κ_SRME.weight = 0.3
    CPS_CONFIG.RMSD.weight = 0.1

    // Find and click the reset button
    const reset_button = document.body.querySelector(`.reset-button`) as HTMLButtonElement
    expect(reset_button).toBeDefined()

    reset_button.click()

    // Check that weights were reset to default values
    expect(CPS_CONFIG.F1.weight).toBe(DEFAULT_CPS_CONFIG.F1.weight)
    expect(CPS_CONFIG.κ_SRME.weight).toBe(DEFAULT_CPS_CONFIG.κ_SRME.weight)
    expect(CPS_CONFIG.RMSD.weight).toBe(DEFAULT_CPS_CONFIG.RMSD.weight)

    // Check that update_models_cps was called
    expect(update_models_cps).toHaveBeenCalled()
  })

  it(`displays correct weight percentages`, () => {
    mount(RadarChart, { target: document.body })

    // Check that the weight percentages are displayed correctly
    // Since we use foreignObject with HTML, we need to look for small elements
    const weight_values = Array.from(
      document.body.querySelectorAll(`svg foreignObject small`),
    ).map((el) => el.textContent)

    // Should match the default weights from DEFAULT_CPS_CONFIG
    const expected_values = Object.values(DEFAULT_CPS_CONFIG).map(
      (part) => `${((part.weight as number) * 100).toFixed(0)}%`,
    )

    expect(weight_values).toEqual(expected_values)
  })

  it(`renders axis labels with correct content`, () => {
    mount(RadarChart, { target: document.body })

    // Check that axis labels are rendered using foreignObject
    const foreign_objects = document.body.querySelectorAll(`svg foreignObject`)
    expect(foreign_objects.length).toBe(3)

    // Get the text content from foreignObject elements
    const label_texts = Array.from(foreign_objects).map((el) => el.textContent?.trim())

    // Verify that each label contains the respective metric name and percentage
    expect(label_texts).toEqual([`F1 50%`, `κSRME 40%`, `RMSD 10%`])
  })

  it(`renders tooltip with config description`, () => {
    mount(RadarChart, {
      target: document.body,
    })

    // Check that the info icon exists (which is part of the tooltip)
    const info_icon = document.body.querySelector(`svg[data-title="Info"]`)
    expect(info_icon).toBeDefined()
  })

  it(`renders correct visualization elements`, () => {
    mount(RadarChart, { target: document.body })

    // Check for the triangle path
    const triangle = document.body.querySelector(`path[d^="M"][d$="Z"]`)
    expect(triangle).toBeDefined()

    // Check for the colored metric areas (should be 3 paths)
    const colored_areas = document.body.querySelectorAll(
      `path[fill^="rgb"][opacity="0.5"]`,
    )
    expect(colored_areas.length).toBe(3)

    // Check for the grid circles
    const grid_circles = document.body.querySelectorAll(`circle[fill="none"]`)
    expect(grid_circles.length).toBe(4) // Should be 4 grid circles
    // Verify stroke properties for each grid circle
    grid_circles.forEach((circle) => {
      expect(circle.getAttribute(`stroke`)).toBe(`var(--border)`)
    })
  })

  it(`displays the metric name`, () => {
    mount(RadarChart, { target: document.body })

    const metric_name = document.body.querySelector(`.metric-name`)
    expect(metric_name).toBeDefined()
    expect(metric_name?.textContent?.trim()).toContain(ALL_METRICS.CPS.key)
  })

  it(`renders reset button with correct attributes`, () => {
    mount(RadarChart, { target: document.body })

    const reset_button = document.body.querySelector(`.reset-button`)
    expect(reset_button).toBeDefined()
    expect(reset_button?.getAttribute(`title`)).toBe(`Reset to default weights`)
    expect(reset_button?.textContent?.trim()).toBe(`Reset`)
  })
})
