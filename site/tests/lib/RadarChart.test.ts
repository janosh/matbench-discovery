import { MODELS } from '$lib'
import RadarChart from '$lib/plot/RadarChart.svelte'
import app_css from '../../src/app.css?raw'
import {
  CPS_CONFIG,
  DEFAULT_CPS_CONFIG,
  CDS_CONFIG,
  DEFAULT_CDS_CONFIG,
} from '$lib/combined-scores.svelte'
import { ALL_METRICS } from '$lib/labels'
import { update_models_cps } from '$lib/models.svelte'
import { format_num } from 'matterviz'
import { flushSync } from 'svelte'
import { describe, expect, it, vi } from 'vitest'
import { doc_query, mount } from '../index'

// Mock the update_models_cps function to avoid side effects
vi.mock(`$lib/models.svelte`, async () => {
  const actual = await vi.importActual(`$lib/models.svelte`)
  return { ...actual, update_models_cps: vi.fn() }
})

describe(`RadarChart`, () => {
  it(`renders with default props`, () => {
    mount(RadarChart, { target: document.body })

    // Check that the component rendered with an SVG element
    const svg = document.querySelector(`svg[aria-label^="Radar chart"]`)
    if (!(svg instanceof SVGSVGElement)) throw new Error(`Radar chart SVG not found`)
    expect(svg.getAttribute(`viewBox`)).toBe(`0 0 200 200`)
    expect(svg.getAttribute(`aria-label`)).toContain(`adjusting metric weights`)

    // Check that all three axis lines are rendered (for F1, kappa, RMSD)
    const axis_lines = svg.querySelectorAll(`line`)
    expect(axis_lines).toHaveLength(3)
    // Verify stroke properties for each axis line
    axis_lines.forEach((line) => {
      expect(line.getAttribute(`stroke`)).toBe(`var(--border)`)
      expect(line.getAttribute(`stroke-width`)).toBe(`1`)
    })

    // Check that the triangle area is rendered
    const triangle_area = svg.querySelector(`path[fill="var(--nav-bg)"]`)
    if (!(triangle_area instanceof SVGPathElement))
      throw new Error(`Triangle path not found`)
    expect(triangle_area.getAttribute(`stroke`)).toBe(`var(--border)`)

    // Check that the draggable point is rendered
    expect(svg.querySelectorAll(`circle[role="button"]`)).toHaveLength(2)
  })

  it(`accepts size prop`, () => {
    const custom_size = 300

    mount(RadarChart, { target: document.body, props: { size: custom_size } })

    const svg = document.querySelector(`svg[aria-label^="Radar chart"]`)
    if (!(svg instanceof SVGSVGElement)) throw new Error(`Radar chart SVG not found`)
    expect(svg.getAttribute(`viewBox`)).toBe(`0 0 ${custom_size} ${custom_size}`)
    expect(svg.getAttribute(`width`)).toBe(String(custom_size))
    expect(svg.getAttribute(`height`)).toBe(String(custom_size))

    // Verify that at least one circle is rendered, which indirectly confirms
    // that the component initialized properly with the size prop
    expect(svg.querySelectorAll(`circle`)).toHaveLength(6)
  })

  it(`resets weights when reset button is clicked`, () => {
    mount(RadarChart, { target: document.body })

    // The reset control stays hidden until the user moves the knob.
    expect(document.querySelector(`.reset-button`)).toBeNull()
    doc_query<SVGSVGElement>(`svg[aria-label^="Radar chart"]`).dispatchEvent(
      new MouseEvent(`click`, { bubbles: true, clientX: 100, clientY: 100 }),
    )
    flushSync()

    const reset_button = doc_query<HTMLButtonElement>(`.reset-button`)
    expect(reset_button.getAttribute(`aria-label`)).toBe(`Reset to default weights`)
    expect(reset_button.querySelector(`svg`)).toBeInstanceOf(SVGSVGElement)
    expect(reset_button.parentElement?.classList.contains(`chart-title`)).toBe(true)

    reset_button.focus()
    reset_button.click()
    flushSync()

    // Check that weights were reset to default values
    expect(CPS_CONFIG.F1.weight).toBe(DEFAULT_CPS_CONFIG.F1.weight)
    expect(CPS_CONFIG.κ_SRME.weight).toBe(DEFAULT_CPS_CONFIG.κ_SRME.weight)
    expect(CPS_CONFIG.RMSD.weight).toBe(DEFAULT_CPS_CONFIG.RMSD.weight)

    // Check that update_models_cps was called
    expect(update_models_cps).toHaveBeenCalledWith(MODELS, CPS_CONFIG)
    expect(document.querySelector(`.reset-button`)).toBeNull()
    expect(document.activeElement).toBe(doc_query(`.metric-name`))
  })

  it(`displays correct weight percentages`, () => {
    mount(RadarChart, { target: document.body })

    // Check that the weight percentages are displayed correctly
    const weight_values = [...document.querySelectorAll(`svg foreignObject`)]
      .map((el) => el.querySelector(`small`)?.textContent)
      .filter(Boolean)

    // mirrors the component's format_num(weight, `.0%`) rendering
    const expected_values = Object.values(DEFAULT_CPS_CONFIG).map((part) =>
      format_num(part.weight, `.0%`),
    )

    expect(weight_values).toStrictEqual(expected_values)
  })

  it(`renders axis labels with correct content`, () => {
    mount(RadarChart, { target: document.body })

    // Check that axis labels are rendered using foreignObject
    const foreign_objects = document.querySelectorAll(`svg foreignObject`)
    expect(foreign_objects).toHaveLength(3)

    // Get the text content from foreignObject elements
    const label_texts = [...foreign_objects].map((el) => el.textContent?.trim())

    // Verify that each label contains the respective metric name and percentage
    expect(label_texts).toStrictEqual([`F1 50%`, `κSRME 40%`, `RMSD 10%`])
  })

  it(`renders correct visualization elements`, () => {
    mount(RadarChart, { target: document.body })

    // Info icon is part of the metric-name tooltip
    expect(document.querySelector(`.metric-name svg`)).toBeInstanceOf(SVGSVGElement)

    // Check for the triangle path
    expect(document.querySelector(`path[d^="M"][d$="Z"]`)).toBeInstanceOf(SVGPathElement)

    // Check for the colored metric areas (should be 3 paths)
    const colored_areas = document.querySelectorAll(`path[fill^="rgb"][opacity="0.5"]`)
    expect(colored_areas).toHaveLength(3)

    // Check for the grid circles
    const grid_circles = document.querySelectorAll(`circle[fill="none"]`)
    expect(grid_circles).toHaveLength(4) // Should be 4 grid circles
    // Verify stroke properties for each grid circle
    grid_circles.forEach((circle) => {
      expect(circle.getAttribute(`stroke`)).toBe(`var(--border)`)
    })
  })

  it(`displays the metric name`, () => {
    mount(RadarChart, { target: document.body })

    const metric_name = doc_query(`.metric-name`)
    expect(metric_name.textContent?.trim()).toContain(ALL_METRICS.CPS.key)
  })

  it(`keeps the knob at the drop point and ignores the trailing click after a drag`, () => {
    mount(RadarChart, { target: document.body })
    flushSync() // let the bind:this effect set svg_element before dragging
    const svg = doc_query<SVGSVGElement>(`svg[aria-label^="Radar chart"]`)
    const knob = doc_query<SVGCircleElement>(`circle[role="button"]`, svg)
    const read = () => ({
      cx: Number(knob.getAttribute(`cx`)),
      cy: Number(knob.getAttribute(`cy`)),
    })

    // press the knob, drag to an interior point, release (size=200 -> center (100,100))
    knob.dispatchEvent(new MouseEvent(`mousedown`, { bubbles: true, cancelable: true }))
    globalThis.dispatchEvent(
      new MouseEvent(`mousemove`, {
        bubbles: true,
        cancelable: true,
        clientX: 90,
        clientY: 110,
      }),
    )
    flushSync()
    const dropped = read()
    expect(dropped.cx).toBeCloseTo(90, 0)
    expect(dropped.cy).toBeCloseTo(110, 0)

    globalThis.dispatchEvent(
      new MouseEvent(`mouseup`, { bubbles: true, cancelable: true }),
    )
    flushSync()
    expect(read()).toEqual(dropped) // drag-end round-trip must not move the knob

    // the browser fires a trailing click after a drag; with wildly different coords
    // (mimicking a layout shift from the weight update) it must NOT move the knob
    svg.dispatchEvent(
      new MouseEvent(`click`, {
        bubbles: true,
        cancelable: true,
        clientX: 160,
        clientY: 100,
      }),
    )
    flushSync()
    expect(read()).toEqual(dropped)

    // a genuine standalone click afterwards should still move the knob
    svg.dispatchEvent(
      new MouseEvent(`click`, {
        bubbles: true,
        cancelable: true,
        clientX: 100,
        clientY: 100,
      }),
    )
    flushSync()
    expect(read()).not.toEqual(dropped)
  })

  it(`supports 4-corner configs: weights stay non-negative and sum to 1 for clicks inside and outside the polygon`, () => {
    mount(RadarChart, {
      target: document.body,
      props: {
        config: CDS_CONFIG,
        default_config: DEFAULT_CDS_CONFIG,
        title_label: ALL_METRICS.diatomics_combined_score,
        on_change: () => {},
      },
    })
    flushSync()
    const svg = doc_query<SVGSVGElement>(`svg[aria-label^="Radar chart"]`)
    const click_at = (x: number, y: number) => {
      svg.dispatchEvent(
        new MouseEvent(`click`, {
          bubbles: true,
          cancelable: true,
          clientX: x,
          clientY: y,
        }),
      )
      flushSync()
      return Object.values(CDS_CONFIG).map(({ weight }) => weight)
    }

    // one axis line + label per pillar
    expect(svg.querySelectorAll(`line`)).toHaveLength(4)
    expect(svg.querySelectorAll(`foreignObject`)).toHaveLength(4)

    // interior click, off-center + outside-polygon click (must clamp to boundary):
    // Wachspress coordinates guarantee weights >= 0 summing to 1 either way
    for (const [x, y] of [
      [130, 80],
      [2, 3],
    ]) {
      const weights = click_at(x, y)
      expect(weights.every((weight) => weight >= 0)).toBe(true)
      expect(weights.reduce((sum, weight) => sum + weight, 0)).toBeCloseTo(1, 10)
    }

    const corner_0_offset = 80 * Math.SQRT1_2

    // corner-adjacent click: the nearest pillar dominates (≈156.6 for size 200)
    const corner = 100 + corner_0_offset
    const [first_pillar, ...rest] = click_at(corner, corner)
    expect(first_pillar).toBeGreaterThan(0.9)
    expect(Math.max(...rest)).toBeLessThan(0.1)

    // the reset/default knob position must round-trip: DEFAULT_CDS_CONFIG weights
    // are canonical (bilinear) coordinates of the point 1/3 of the way to the
    // accuracy corner, so clicking exactly there must reproduce the defaults -
    // guards against defaults the knob can't express (which would make the reset
    // knob position visually disagree with the displayed percentages). Last click,
    // so it doubles as the shared-module-state restore for other tests
    const defaults = Object.values(DEFAULT_CDS_CONFIG).map(({ weight }) => weight)
    const reset_pos = 100 + corner_0_offset / 3
    for (const [idx, weight] of click_at(reset_pos, reset_pos).entries()) {
      expect(weight).toBeCloseTo(defaults[idx], 6)
    }
  })

  it(`shows the reset button for non-default weights even when the knob lands on the default position`, () => {
    // 50/0/50/0 on opposite corners of a 4-corner chart has its weighted centroid at
    // the polygon center - for these configs the weights->point map is many-to-one,
    // so reset-button visibility must derive from the weights, not knob geometry
    for (const [key, pillar] of Object.entries(CDS_CONFIG)) {
      pillar.weight = [`accuracy`, `speed`].includes(key) ? 0.5 : 0
    }
    try {
      mount(RadarChart, {
        target: document.body,
        props: {
          config: CDS_CONFIG,
          default_config: DEFAULT_CDS_CONFIG,
          title_label: ALL_METRICS.diatomics_combined_score,
          on_change: () => {},
        },
      })
      flushSync()

      doc_query<HTMLButtonElement>(`.reset-button`).click()
      flushSync()

      expect(Object.values(CDS_CONFIG).map(({ weight }) => weight)).toEqual(
        Object.values(DEFAULT_CDS_CONFIG).map(({ weight }) => weight),
      )
      expect(document.querySelector(`.reset-button`)).toBeNull()
    } finally {
      // restore shared module state for other tests even if assertions throw
      for (const [key, pillar] of Object.entries(DEFAULT_CDS_CONFIG)) {
        CDS_CONFIG[key as keyof typeof CDS_CONFIG].weight = pillar.weight
      }
    }
  })

  it(`disables scroll anchoring so adjusting CPS weights can't jump the page`, () => {
    // Dragging the knob re-renders the metrics table, reflowing content above the
    // viewport. Without overflow-anchor:none the browser scrolls to keep an anchor
    // element in place (page jumps to the scatter plot). happy-dom can't exercise
    // scroll anchoring, so guard the CSS rule that fixes it directly.
    expect(app_css).toMatch(/html\s*\{[^}]*overflow-anchor:\s*none/)
  })
})
