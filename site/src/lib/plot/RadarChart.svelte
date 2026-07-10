<script lang="ts">
  import type { CpsConfig } from '$lib/combined-scores.svelte'
  import { ALL_METRICS } from '$lib/labels'
  import type { Label } from '$lib/types'
  import { format_num, Icon, type Point } from 'matterviz'
  import { tooltip } from 'svelte-multiselect/attachments'
  import { CPS_CONFIG, DEFAULT_CPS_CONFIG } from '$lib/combined-scores.svelte'
  import { MODELS, update_models_cps } from '$lib/models.svelte'

  // any weighted score with >= 3 components works (CPS is the default; CMDS and CDS
  // use the same UI with 3 and 4 corners respectively)
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
  // cursor against a moved rect and jump the knob (often into a polygon corner).
  // Skip that one click.
  let suppress_next_click = false
  let svg_element = $state<SVGSVGElement | null>(null)
  let title_element = $state<HTMLElement | null>(null)
  let radius = $derived(size / 2)
  let center = $derived({ x: radius, y: radius })

  const colors = [
    `rgb(255, 99, 132)`, // red
    `rgb(255, 206, 86)`, // yellow
    `rgb(54, 162, 235)`, // blue
    `rgb(75, 192, 120)`, // green
    `rgb(153, 102, 255)`, // purple
    `rgb(255, 159, 64)`, // orange
  ]

  // Corner angle for axis idx: when the corner count is divisible by 4, rotate by
  // half a sector so no corner points straight up into the title text
  function corner_angle(idx: number, n_corners: number): number {
    const offset = n_corners % 4 === 0 ? Math.PI / n_corners : 0
    return (2 * Math.PI * idx) / n_corners + offset
  }

  // Compute axes points coordinates
  let axis_points = $derived(
    Object.values(config).map((_, idx, { length }) => {
      const angle = corner_angle(idx, length)
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
  // svelte-ignore state_referenced_locally
  let point = $state<Point>(point_from_weights(config))
  // Compare weights, not knob geometry: for 4+ corners the weights->point map is
  // many-to-one (e.g. 50/0/50/0 on opposite corners lands on the center, same as
  // equal weights), so a position check would hide the Reset button for non-default
  // weights. Weight-based gating also keeps the button stable mid-drag (weights only
  // update on drop). Keys absent from default_config don't count as custom since
  // reset_weights can't restore them anyway.
  let has_custom_weights = $derived(
    Object.keys(config).some((key) => {
      const default_weight = default_config[key]?.weight
      return (
        default_weight !== undefined &&
        Math.abs(config[key].weight - default_weight) > 1e-6
      )
    }),
  )

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
    title_element?.focus({ preventScroll: true })
  }

  // Derive weights from the knob position via Wachspress generalized barycentric
  // coordinates: non-negative inside any convex polygon, sum to 1 by construction
  // (so the knob can never express weights totaling more or less than 1), reduce to
  // classic barycentric coordinates for triangles, and satisfy Σ λᵢ·vᵢ = p so they
  // round-trip exactly with point_from_weights. Every point inside the polygon is a
  // valid weight assignment - the reachable set of weighted centroids IS the
  // straight-edged polygon, no curved boundary needed. For n > 3 corners the map
  // point -> weights is no longer unique (a 2D point can't encode an
  // (n-1)-simplex); Wachspress picks the canonical smooth choice.
  function update_weights_from_point() {
    const keys = Object.keys(config)
    const n_corners = axis_points.length
    const wachspress = (pt: Point) =>
      axis_points.map((vertex, idx) => {
        const prev = axis_points[(idx - 1 + n_corners) % n_corners]
        const next = axis_points[(idx + 1) % n_corners]
        const corner_area = calc_triangle_area(prev, vertex, next)
        return (
          corner_area /
          (calc_triangle_area(pt, prev, vertex) * calc_triangle_area(pt, vertex, next))
        )
      })
    let raw = wachspress(point)
    // a knob exactly on a corner/edge zeroes Wachspress denominators (-> Infinity):
    // recompute nudged sub-pixel toward the polygon center where the limit weights
    // emerge numerically. Interior points keep the exact round-trip.
    if (!raw.every(Number.isFinite)) {
      raw = wachspress({
        x: point.x + (center.x - point.x) * 1e-6,
        y: point.y + (center.y - point.y) * 1e-6,
      })
    }
    const total = raw.reduce((sum, val) => sum + val, 0)
    let new_values = raw.map((val) => val / total)

    // Snap to equal weights when very close to the center
    const dist_from_center = Math.hypot(point.x - center.x, point.y - center.y)
    if (dist_from_center < radius * 0.05) new_values = keys.map(() => 1 / n_corners)

    // Update weights in the config directly (axis order == config key order)
    for (const [idx, key] of keys.entries()) config[key].weight = new_values[idx]
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
  function clamp_to_polygon(x: number, y: number): Point {
    const pt = { x, y }
    return is_point_in_polygon(pt) ? pt : closest_point_on_polygon(pt)
  }

  // Move the point to a position, constraining to the corner polygon if needed
  function move_to_position(click_x: number, click_y: number) {
    point = clamp_to_polygon(click_x, click_y)
    update_weights_from_point()
  }

  // Visually move the knob during a drag without touching weights (see end_drag)
  function move_point_to_position(x: number, y: number) {
    point = clamp_to_polygon(x, y)
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
    if (!is_dragging) return
    is_dragging = false
    // Update weights only on drag end: doing it live would rerender the table
    // and scroll the viewport mid-drag
    update_weights_from_point()
    globalThis.removeEventListener(`mousemove`, handle_drag)
    globalThis.removeEventListener(`touchmove`, handle_drag)
    globalThis.removeEventListener(`mouseup`, end_drag)
    globalThis.removeEventListener(`touchend`, end_drag)
  }

  // Convex polygon with consistently ordered vertices: pt is inside iff it lies on
  // the same side of every edge
  function is_point_in_polygon(pt: Point): boolean {
    let sign = 0
    for (const [idx, v1] of axis_points.entries()) {
      const v2 = axis_points[(idx + 1) % axis_points.length]
      const cross = (v2.x - v1.x) * (pt.y - v1.y) - (v2.y - v1.y) * (pt.x - v1.x)
      if (cross === 0) continue
      if (sign === 0) sign = Math.sign(cross)
      else if (Math.sign(cross) !== sign) return false
    }
    return true
  }

  // Nearest boundary point: the closest of the per-edge segment projections
  function closest_point_on_polygon(pt: Point): Point {
    let best = axis_points[0]
    let best_dist = Infinity
    for (const [idx, v1] of axis_points.entries()) {
      const v2 = axis_points[(idx + 1) % axis_points.length]
      const candidate = closest_point_on_line(pt, v1, v2)
      const dist = (candidate.x - pt.x) ** 2 + (candidate.y - pt.y) ** 2
      if (dist < best_dist) [best_dist, best] = [dist, candidate]
    }
    return best
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

<div class="radar-chart" style="--chart-size: {size}px">
  <div class="chart-title">
    <span
      bind:this={title_element}
      class="metric-name"
      title={title_label.description}
      tabindex="-1"
      {@attach tooltip()}
    >
      {@html title_label.label ?? title_label.key}
      <Icon icon="Info" />
    </span>
    {#if has_custom_weights}
      <button
        class="reset-button"
        onclick={reset_weights}
        aria-label="Reset to default weights"
      >
        <Icon icon="Reset" />
      </button>
    {/if}
  </div>

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
      {@const angle = corner_angle(idx, Object.values(config).length)}
      {@const { x, y } = axis_points[idx]}
      <!-- push labels above the center further out so they don't collide with the
      title (SVG y grows downward, so sin < 0 = upper half) -->
      {@const label_radius = radius * (Math.sin(angle) < -1e-9 ? 1.0 : 0.9)}
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
        style:color={colors[idx % colors.length]}
      >
        {@html weight.label ?? weight.key}
        <small>{format_num(weight.weight, `.0%`)}</small>
      </foreignObject>
    {/each}

    <!-- Draggable corner-polygon area -->
    <path
      d="{axis_points
        .map(({ x, y }, idx) => `${idx === 0 ? `M` : `L`} ${x} ${y}`)
        .join(` `)} Z"
      fill="var(--nav-bg)"
      stroke="var(--border)"
      stroke-width="1"
    />

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
    {#each axis_points as { x: x0, y: y0 }, idx (idx)}
      {@const path = `M ${center.x} ${center.y} L ${x0} ${y0} L ${point.x} ${point.y} Z`}
      <path d={path} fill={colors[idx % colors.length]} stroke="none" opacity="0.5" />
    {/each}

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
    display: grid;
    grid-template-columns: 3ex var(--chart-size) 3ex;
    justify-content: center;
    justify-items: center;
    width: calc(var(--chart-size) + 6ex);
    margin: 0 auto;
    background: var(--card-bg);
    border-radius: 4px;
  }
  .chart-title,
  svg {
    grid-column: 2;
  }
  svg {
    touch-action: none; /* Prevents default touch behaviors */
    cursor: pointer; /* Show pointer cursor to hint clickability */
    overflow: visible;
  }
  svg:focus:not(:focus-visible) {
    outline: none;
  }
  svg:focus-visible {
    outline: 2px solid var(--link-color);
    outline-offset: 2px;
  }
  .chart-title {
    display: flex;
    align-items: center;
    gap: 0.5em;
    min-height: 1em;
    white-space: nowrap;
  }
  .metric-name {
    display: inline-flex;
    align-items: center;
    gap: 0.4em;
    line-height: 1;
  }
  .metric-name:focus {
    outline: 2px solid var(--link-color);
    outline-offset: 2px;
  }
  .reset-button {
    display: inline-flex;
    align-items: center;
    background: transparent;
    color: inherit;
    padding: 0.2em;
  }
</style>
