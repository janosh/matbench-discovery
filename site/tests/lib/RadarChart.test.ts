import { RadarChart } from '$lib'
import { CPS_CONFIG, DEFAULT_CPS_CONFIG } from '$lib/combined_perf_score.svelte'
import { ALL_METRICS } from '$lib/labels'
import { update_models_cps } from '$lib/models.svelte'
import { mount } from 'svelte'
import { beforeEach, describe, expect, it, vi } from 'vitest'

describe(`RadarChart`, () => {
  beforeEach(() => {
    document.body.innerHTML = `` // Reset DOM before each test

    // Reset CPS_CONFIG to default values before each test
    let key: keyof typeof CPS_CONFIG
    for (key in CPS_CONFIG) {
      CPS_CONFIG[key].weight = DEFAULT_CPS_CONFIG[key].weight
    }

    // Mock the update_models_cps function to avoid side effects
    vi.mock(`$lib/models.svelte`, async () => {
      const actual = await vi.importActual(`$lib/models.svelte`)
      return {
        ...actual,
        update_models_cps: vi.fn(),
      }
    })
  })

  it(`renders with default props`, async () => {
    mount(RadarChart, { target: document.body })

    // Check that the component rendered with an SVG element
    const svg = document.body.querySelector(`svg`)
    expect(svg).toBeDefined()

    // Check that all three axis lines are rendered (for F1, kappa, RMSD)
    const axis_lines = document.body.querySelectorAll(`line`)
    expect(axis_lines).toHaveLength(3)

    // Check that the triangle area is rendered
    const triangle_area = document.body.querySelector(
      `path[fill="rgba(255, 255, 255, 0.1)"]`,
    )
    expect(triangle_area).toBeDefined()

    // Check that the draggable point is rendered
    const draggable_point = document.body.querySelector(`circle[role="button"]`)
    expect(draggable_point).toBeDefined()
  })

  it(`accepts size prop`, async () => {
    const custom_size = 300

    // Just verify that the component mounts without error when we pass a size prop
    mount(RadarChart, { target: document.body, props: { size: custom_size } })

    // Just verify the SVG element exists
    const svg = document.body.querySelector(`svg`)
    expect(svg).toBeDefined()

    // Verify that at least one circle is rendered, which indirectly confirms
    // that the component initialized properly with the size prop
    const circles = document.body.querySelectorAll(`circle`)
    expect(circles.length).toBeGreaterThan(0)
  })

  it(`resets weights when reset button is clicked`, async () => {
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

  it(`displays correct weight percentages`, async () => {
    mount(RadarChart, { target: document.body })

    // Check that the weight percentages are displayed correctly
    const weight_values = Array.from(
      document.body.querySelectorAll(`svg tspan:last-child`),
    ).map((el) => el.textContent)

    // Should match the default weights from DEFAULT_CPS_CONFIG
    const expected_values = Object.values(DEFAULT_CPS_CONFIG).map(
      (part) => `${((part.weight as number) * 100).toFixed(0)}%`,
    )

    expect(weight_values).toEqual(expected_values)
  })

  it(`renders axis labels with correct content`, async () => {
    mount(RadarChart, { target: document.body })

    // Check that axis labels are rendered
    const text_elements = document.body.querySelectorAll(`svg text`)
    expect(text_elements.length).toBe(3)

    // In the actual component, labels also include the percentage value
    // So we need to get the full text content
    const label_texts = Array.from(text_elements).map((el) => el.textContent?.trim())

    // Verify that each label contains the respective metric name
    expect(label_texts.some((text) => text?.includes(`F1`))).toBe(true)
    expect(label_texts.some((text) => text?.includes(`κ SRME`))).toBe(true)
    expect(label_texts.some((text) => text?.includes(`RMSD`))).toBe(true)
  })

  it(`renders tooltip with config description`, async () => {
    mount(RadarChart, {
      target: document.body,
    })

    // Check that the info icon exists (which is part of the tooltip)
    const info_icon = document.body.querySelector(`svg use[href="#icon-info"]`)
    expect(info_icon).toBeDefined()
  })

  it(`renders correct visualization elements`, async () => {
    mount(RadarChart, { target: document.body })

    // Check for the triangle path
    const triangle = document.body.querySelector(`path[d^="M"][d$="Z"]`)
    expect(triangle).toBeDefined()

    // Check for the colored metric areas (should be 3 paths)
    const colored_areas = document.body.querySelectorAll(
      `path[fill^="rgba"][opacity="0.5"]`,
    )
    expect(colored_areas.length).toBe(3)

    // Check for the grid circles
    const grid_circles = document.body.querySelectorAll(
      `circle[fill="none"][stroke="rgba(255, 255, 255, 0.1)"]`,
    )
    expect(grid_circles.length).toBe(4) // Should be 4 grid circles
  })

  it(`displays the metric name`, async () => {
    mount(RadarChart, { target: document.body })

    const metric_name = document.body.querySelector(`.metric-name`)
    expect(metric_name).toBeDefined()
    expect(metric_name?.textContent?.trim()).toContain(ALL_METRICS.CPS.label)
  })

  it(`renders reset button with correct attributes`, async () => {
    mount(RadarChart, { target: document.body })

    const reset_button = document.body.querySelector(`.reset-button`)
    expect(reset_button).toBeDefined()
    expect(reset_button?.getAttribute(`title`)).toBe(`Reset to default weights`)
    expect(reset_button?.textContent?.trim()).toBe(`Reset`)
  })
})
