<script lang="ts">
  import type { CpsConfig } from '$lib/combined_perf_score.svelte'
  import { ALL_METRICS } from '$lib/labels'
  import type { Label } from '$lib/types'
  import { format_num, Icon, type Point } from 'matterviz'
  import { tooltip } from 'svelte-multiselect/attachments'
  import { CPS_CONFIG, DEFAULT_CPS_CONFIG } from '$lib/combined_perf_score.svelte'
  import { MODELS, update_models_cps } from '$lib/models.svelte'

  // any 3-component weighted score works (CPS is the default, CMDS uses the same UI)
  type WeightsConfig = Record<string, Label & { weight: number }>

  let {
    size = 200,
    config = CPS_CONFIG,
    default_config = DEFAULT_CPS_CONFIG,
    title_label = ALL_METRICS.CPS,
    on_change = (cfg: WeightsConfig) => update_models_cps(MODELS, cfg as CpsConfig),
  }: {
    size?: number
    config?: WeightsConfig
    default_config?: WeightsConfig
    title_label?: Label
    on_change?: (config: WeightsConfig) => void
  } = $props()

  // State for the draggable point
  let is_dragging = $state(false)
  // A pointer interaction that starts on the knob makes the browser fire a
  // trailing click on the SVG. Dropping the knob updates the weights and can
  // reflow the page (shifting the SVG), so handling that click would remap the
  // cursor against a moved rect and jump the knob (often into a triangle
  // corner). Skip that one click.
  let suppress_next_click = false
  let svg_element = $state<SVGSVGElement | null>(null)
  let radius = $derived(size / 2)
  let center = $derived({ x: radius, y: radius })

  const colors = [
    `rgb(255, 99, 132)`, // Red for F1
    `rgb(255, 206, 86)`, // Yellow for kappa
    `rgb(54, 162, 235)`, // Blue for RMSD
  ]

  // Compute axes points coordinates
  let axis_points = $derived(
    Object.values(config).map((_, idx) => {
      const angle = (2 * Math.PI * idx) / Object.values(config).length
      const x = center.x + Math.cos(angle) * radius * 0.8
      const y = center.y + Math.sin(angle) * radius * 0.8
      return { x, y }
    }),
  )

  // Knob position from weights via barycentric coordinates (pure math, no DOM)
  function point_from_weights(weights: WeightsConfig): Point {
    let { x, y } = center
    Object.values(weights).forEach(({ weight }, idx) => {
      x += (axis_points[idx].x - center.x) * weight
      y += (axis_points[idx].y - center.y) * weight
    })
    return { x, y }
  }

  // Initialized eagerly (not in an effect) so the knob renders at the correct
  // position on first paint and in prerendered HTML instead of jumping from (0, 0)
  let point = $state<Point>(point_from_weights(config))

  // Keep knob + model scores in sync when weights change (drag end, reset,
  // or edits from elsewhere)
  $effect(() => {
    point = point_from_weights(config)
    on_change(config)
  })

  // Reset to initial weights (the effect above re-derives knob position + scores).
  // Skip keys absent from default_config so a divergent config pair can't crash.
  function reset_weights() {
    for (const key of Object.keys(config) as (keyof typeof config)[]) {
      const default_weight = default_config[key]?.weight
      if (default_weight !== undefined) config[key].weight = default_weight
    }
  }

  // Derive weights from the knob position via barycentric coordinates
  function update_weights_from_point() {
    if (Object.values(config).length !== 3) {
      console.error(`This implementation only supports exactly 3 metrics/dimensions`)
      return
    }

    const [a, b, c] = axis_points
    const triangle_area = calc_triangle_area(a, b, c)
    const area1 = calc_triangle_area(point, b, c)
    const area2 = calc_triangle_area(point, a, c)
    const area3 = calc_triangle_area(point, a, b)

    let new_values = [area1 / triangle_area, area2 / triangle_area, area3 / triangle_area]

    // Snap to equal weights when very close to the center
    const dist_from_center = Math.hypot(point.x - center.x, point.y - center.y)
    if (dist_from_center < radius * 0.05) {
      new_values = [1 / 3, 1 / 3, 1 / 3]
    }

    // Update weights in the config directly (axis order == config key order)
    for (const [idx, key] of Object.keys(config).entries()) {
      config[key].weight = new_values[idx]
    }
  }

  // Helper to calculate triangle area using cross product
  const calc_triangle_area = (p1: Point, p2: Point, p3: Point) =>
    Math.abs((p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y)) / 2)

  // Handle click on SVG to jump to position
  function handle_svg_click(event: MouseEvent) {
    event.preventDefault()

    // Ignore the synthetic click the browser fires right after a knob drag
    if (suppress_next_click) {
      suppress_next_click = false
      return
    }

    const target = event.currentTarget as SVGSVGElement
    const svg_rect = target.getBoundingClientRect()
    move_to_position(event.clientX - svg_rect.left, event.clientY - svg_rect.top)
  }

  // Handle keyboard activation - move to SVG center (equal weights)
  function handle_keyboard_click(event: KeyboardEvent) {
    if (event.key === `Enter` || event.key === ` `) {
      event.preventDefault()
      move_to_position(center.x, center.y)
    }
  }

  // Shared by the click and drag paths so their clamping behavior can't drift
  function clamp_to_triangle(x: number, y: number): Point {
    const [a, b, c] = axis_points
    const pt = { x, y }
    return Object.values(config).length === 3 && is_point_in_triangle(pt, a, b, c)
      ? pt
      : get_closest_point_on_triangle(pt, a, b, c)
  }

  // Move the point to a position, constraining to triangle if needed
  function move_to_position(click_x: number, click_y: number) {
    point = clamp_to_triangle(click_x, click_y)
    update_weights_from_point()
  }

  // Visually move the knob during a drag without touching weights (see end_drag)
  function move_point_to_position(x: number, y: number) {
    point = clamp_to_triangle(x, y)
  }

  function start_drag(event: MouseEvent | TouchEvent) {
    event.preventDefault()
    is_dragging = true
    suppress_next_click = true

    globalThis.addEventListener(`mousemove`, handle_drag)
    globalThis.addEventListener(`touchmove`, handle_drag, { passive: false })
    globalThis.addEventListener(`mouseup`, end_drag)
    globalThis.addEventListener(`touchend`, end_drag)
  }

  function handle_drag(event: MouseEvent | TouchEvent) {
    if (!is_dragging) return
    event.preventDefault()
    if (!svg_element) return

    const rect = svg_element.getBoundingClientRect()
    const { clientX, clientY } = event instanceof MouseEvent ? event : event.touches[0]
    move_point_to_position(clientX - rect.left, clientY - rect.top)
  }

  function end_drag() {
    if (is_dragging) {
      is_dragging = false
      // Update weights only on drag end: doing it live would rerender the table
      // and scroll the viewport mid-drag
      update_weights_from_point()
      globalThis.removeEventListener(`mousemove`, handle_drag)
      globalThis.removeEventListener(`touchmove`, handle_drag)
      globalThis.removeEventListener(`mouseup`, end_drag)
      globalThis.removeEventListener(`touchend`, end_drag)
    }
  }

  // Helper to check if a point pt is inside a triangle with corners c1, c2, c3
  function is_point_in_triangle(pt: Point, c1: Point, c2: Point, c3: Point) {
    // Compute barycentric coordinates
    const denominator = (c2.y - c3.y) * (c1.x - c3.x) + (c3.x - c2.x) * (c1.y - c3.y)
    const alpha =
      ((c2.y - c3.y) * (pt.x - c3.x) + (c3.x - c2.x) * (pt.y - c3.y)) / denominator
    const beta =
      ((c3.y - c1.y) * (pt.x - c3.x) + (c1.x - c3.x) * (pt.y - c3.y)) / denominator
    const gamma = 1 - alpha - beta

    // If all coordinates are between 0 and 1, point is inside
    return alpha >= 0 && beta >= 0 && gamma >= 0 && alpha <= 1 && beta <= 1 && gamma <= 1
  }

  // Helper to find closest point pt on triangle with corners c1, c2, c3
  function get_closest_point_on_triangle(pt: Point, c1: Point, c2: Point, c3: Point) {
    // First, find the closest point on each edge
    const p12 = closest_point_on_line(pt, c1, c2)
    const p23 = closest_point_on_line(pt, c2, c3)
    const p31 = closest_point_on_line(pt, c3, c1)

    // Then, find which of those points is closest to p
    const dist_p12 = (p12.x - pt.x) ** 2 + (p12.y - pt.y) ** 2
    const dist_p23 = (p23.x - pt.x) ** 2 + (p23.y - pt.y) ** 2
    const dist_p31 = (p31.x - pt.x) ** 2 + (p31.y - pt.y) ** 2

    if (dist_p12 <= dist_p23 && dist_p12 <= dist_p31) return p12
    if (dist_p23 <= dist_p12 && dist_p23 <= dist_p31) return p23
    return p31
  }

  // Helper to find closest point on a line segment
  function closest_point_on_line(p: Point, a: Point, b: Point) {
    const atob = { x: b.x - a.x, y: b.y - a.y }
    const atop = { x: p.x - a.x, y: p.y - a.y }
    const len = atob.x * atob.x + atob.y * atob.y
    let dot = atop.x * atob.x + atop.y * atob.y
    const t = Math.max(0, Math.min(1, dot / len))

    return { x: a.x + atob.x * t, y: a.y + atob.y * t }
  }
</script>

<div class="radar-chart">
  <span class="metric-name" title={title_label.description} {@attach tooltip()}>
    {@html title_label.label ?? title_label.key}
    <Icon icon="Info" />
  </span>

  <button class="reset-button" onclick={reset_weights} title="Reset to default weights">
    Reset
  </button>

  <svg
    bind:this={svg_element}
    width={size}
    height={size}
    viewBox="0 0 {size} {size}"
    onclick={handle_svg_click}
    onkeydown={handle_keyboard_click}
    tabindex="0"
    role="button"
    aria-label="Radar chart for adjusting metric weights. Click to set custom weights. Press Enter or Space to reset to equal weights."
  >
    <!-- Axes -->
    {#each Object.values(config) as weight, idx (weight.label)}
      {@const angle = (2 * Math.PI * idx) / Object.values(config).length}
      {@const x = center.x + Math.cos(angle) * radius * 0.8}
      {@const y = center.y + Math.sin(angle) * radius * 0.8}
      {@const label_radius = idx === 2 ? radius * 1.0 : radius * 0.9}
      {@const label_x = center.x + Math.cos(angle) * label_radius}
      {@const label_y = center.y + Math.sin(angle) * label_radius}
      <line
        x1={center.x}
        y1={center.y}
        x2={x}
        y2={y}
        stroke="var(--border)"
        stroke-width="1"
      />

      <!-- Axis labels -->
      <foreignObject
        x={label_x}
        y={label_y}
        style="font-size: 14px; transform: translate(-1ex, -1em); overflow: visible; white-space: nowrap"
        style:color={colors[idx]}
      >
        {@html weight.label ?? weight.key}
        <small>{format_num(weight.weight, `.0%`)}</small>
      </foreignObject>
    {/each}

    <!-- Triangle area -->
    {#if Object.values(config).length === 3}
      <path
        d="M {axis_points[0].x} {axis_points[0].y} L {axis_points[1].x} {axis_points[1]
          .y} L {axis_points[2].x} {axis_points[2].y} Z"
        fill="var(--nav-bg)"
        stroke="var(--border)"
        stroke-width="1"
      />
    {/if}

    <!-- Background circular grid -->
    {#each [0.2, 0.4, 0.6, 0.8] as grid_radius (grid_radius)}
      <circle
        cx={center.x}
        cy={center.y}
        r={radius * grid_radius}
        fill="none"
        stroke="var(--border)"
        stroke-width="1"
      />
    {/each}

    <!-- Colored areas for each metric -->
    {#if Object.values(config).length === 3}
      {@const [{ x, y }, { x: px, y: py }] = [center, point]}
      {#each axis_points as { x: x0, y: y0 }, idx (idx)}
        {@const path = `M ${x} ${y} L ${x0} ${y0} L ${px} ${py} Z`}
        <path d={path} fill={colors[idx]} stroke="none" opacity="0.5" />
      {/each}
    {/if}

    <!-- Draggable knob: first element is larger invisible hit area for the smaller visible knob above it -->
    <!-- page-bg (not card-bg): the knob needs an opaque fill, card-bg is 0.3-0.4 alpha -->
    {#each [{ fill: `transparent`, r: 20 }, { fill: `var(--page-bg)`, stroke: `var(--text-color)`, r: 7 }] as knob_style (knob_style.r)}
      <circle
        cx={point.x}
        cy={point.y}
        style="cursor: move"
        onmousedown={start_drag}
        ontouchstart={start_drag}
        role="button"
        tabindex="0"
        aria-label="Drag to adjust weight balance"
        {...knob_style}
      />
    {/each}
  </svg>
</div>

<style>
  .radar-chart {
    padding: 1em 3ex 0 0;
    margin: 0;
    position: relative;
    background: var(--card-bg);
    border-radius: 4px;
  }
  svg {
    touch-action: none; /* Prevents default touch behaviors */
    cursor: pointer; /* Show pointer cursor to hint clickability */
    overflow: visible;
  }
  span.metric-name {
    position: absolute;
    top: 2pt;
    left: 40%;
  }
  .reset-button {
    position: absolute;
    top: 4pt;
    right: 5pt;
    background: transparent;
    border: 1px solid var(--border);
    border-radius: 3px;
    cursor: pointer;
    padding: 0.15em 0.35em;
    font-size: 0.8em;
    margin-left: auto;
  }
  .reset-button:hover {
    background: var(--nav-bg);
  }
</style>
