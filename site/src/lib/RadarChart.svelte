<script lang="ts">
  import type { CpsConfig } from '$lib/combined_perf_score.svelte'
  import { METRICS } from '$lib/labels'
  import type { Point } from 'elementari'
  import { Tooltip } from 'svelte-zoo'
  import { CPS_CONFIG, DEFAULT_CPS_CONFIG } from './combined_perf_score.svelte'
  import { update_models_cps } from './models.svelte'

  // Define props interface
  interface Props {
    size?: number
  }
  let { size = 200 }: Props = $props()

  // State for the draggable point
  let is_dragging = $state(false)
  let point = $state<Point>({ x: 0, y: 0 })
  let svg_element: SVGSVGElement
  let radius = size / 2
  let center = { x: radius, y: radius }

  // Reset to initial weights
  function reset_weights() {
    let key: keyof typeof CPS_CONFIG
    for (key in CPS_CONFIG) {
      CPS_CONFIG[key].weight = DEFAULT_CPS_CONFIG[key].weight
    }
    update_point_from_weights(CPS_CONFIG)
    update_models_cps() // Update all model CPS with new weights
  }

  const colors = [
    `rgba(255, 99, 132, 0.7)`, // red for F1
    `rgba(255, 206, 86, 0.7)`, // yellow for kappa
    `rgba(54, 162, 235, 0.7)`, // blue for RMSD
  ]

  // Compute axes points coordinates
  let axis_points = $derived(
    Object.values(CPS_CONFIG).map((_, idx) => {
      const angle = (2 * Math.PI * idx) / Object.values(CPS_CONFIG).length
      const x = center.x + Math.cos(angle) * radius * 0.8
      const y = center.y + Math.sin(angle) * radius * 0.8
      return { x, y }
    }),
  )

  // Initialize point position from weights
  $effect(() => {
    update_point_from_weights(CPS_CONFIG)
    update_models_cps()
  })

  function update_point_from_weights(current_weights: CpsConfig[`parts`]) {
    if (!current_weights || Object.values(current_weights).length < 3) return

    // For 3 axes, we can use barycentric coordinates

    // Calculate weighted position
    let { x, y } = center
    for (let idx = 0; idx < Object.values(current_weights).length; idx++) {
      const weight_value = Object.values(current_weights)[idx].weight as number
      x += (axis_points[idx].x - center.x) * weight_value
      y += (axis_points[idx].y - center.y) * weight_value
    }
    // Update the point with new coordinates
    point = { x, y }
  }

  function update_weights_from_point() {
    // Calculate weights using barycentric coordinates for triangular space
    if (Object.values(CPS_CONFIG).length !== 3) {
      console.error(`This implementation only supports exactly 3 metrics`)
      return
    }

    // Implementation for a 3-point spider diagram
    // Calculate barycentric coordinates of the point relative to the triangle

    // First convert to triangle coordinates
    const [a, b, c] = axis_points
    const triangle_area = calc_triangle_area(a, b, c)
    // Calculate areas of sub-triangles
    const area1 = calc_triangle_area(point, b, c)
    const area2 = calc_triangle_area(point, a, c)
    const area3 = calc_triangle_area(point, a, b)

    // Calculate normalized barycentric weights
    let new_values = [area1 / triangle_area, area2 / triangle_area, area3 / triangle_area]

    // Handle center point specially - equal weights
    const dist_from_center = Math.sqrt(
      Math.pow(point.x - center.x, 2) + Math.pow(point.y - center.y, 2),
    )

    // If very close to center, use equal weights
    if (dist_from_center < radius * 0.05) {
      new_values = [1 / 3, 1 / 3, 1 / 3]
    }

    // Update weights in the CPS_CONFIG directly
    CPS_CONFIG.F1.weight = new_values[0]
    CPS_CONFIG.κ_SRME.weight = new_values[1]
    CPS_CONFIG.RMSD.weight = new_values[2]
  }

  // Helper to calculate triangle area using cross product
  function calc_triangle_area(p1: Point, p2: Point, p3: Point) {
    return Math.abs(
      (p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y)) / 2,
    )
  }

  // Handle click on SVG to jump to position
  function handle_svg_click(event: MouseEvent) {
    // Prevent default behavior including scrolling
    event.preventDefault()

    // Get click coordinates relative to SVG
    const target = event.currentTarget as SVGSVGElement
    const svg_rect = target.getBoundingClientRect()
    const click_x = event.clientX - svg_rect.left
    const click_y = event.clientY - svg_rect.top

    // Check if click is inside the triangle
    const [a, b, c] = axis_points
    const click_point = { x: click_x, y: click_y }

    if (is_point_in_triangle(click_point, a, b, c)) {
      // Update the draggable point position
      point = click_point
    } else {
      // If outside the triangle, find the closest point on the triangle
      point = get_closest_point_on_triangle(click_point, a, b, c)
    }

    // Update weights based on new position
    update_weights_from_point()
  }

  // Move point to a position with triangle constraints
  // during dragging, don't update weights - only move the point visually
  // this prevents table rerendering during drag which causes the viewport to scroll (terrible UX)
  function move_point_to_position(x: number, y: number) {
    const [a, b, c] = axis_points
    if (
      Object.values(CPS_CONFIG).length === 3 &&
      is_point_in_triangle({ x, y }, a, b, c)
    ) {
      point = { x, y }
    } else {
      // If outside the triangle, constrain to the closest point on the triangle
      const closest_point = get_closest_point_on_triangle({ x, y }, a, b, c)
      point = closest_point
    }
  }

  // Handle dragging
  function start_drag(event: MouseEvent | TouchEvent) {
    event.preventDefault()
    is_dragging = true

    // Add global event listeners
    window.addEventListener(`mousemove`, handle_drag)
    window.addEventListener(`touchmove`, handle_drag, { passive: false })
    window.addEventListener(`mouseup`, end_drag)
    window.addEventListener(`touchend`, end_drag)
  }

  function handle_drag(event: MouseEvent | TouchEvent) {
    if (!is_dragging) return
    event.preventDefault()

    // Get SVG coordinates
    const rect = svg_element.getBoundingClientRect()

    // Get position depending on event type
    let client_x: number, client_y: number

    if (event instanceof MouseEvent) {
      client_x = event.clientX
      client_y = event.clientY
    } else {
      client_x = event.touches[0].clientX
      client_y = event.touches[0].clientY
    }

    const x = client_x - rect.left
    const y = client_y - rect.top

    // update point position during drag
    move_point_to_position(x, y)
  }

  function end_drag() {
    if (is_dragging) {
      is_dragging = false
      // Update weights when drag ends
      update_weights_from_point()
      window.removeEventListener(`mousemove`, handle_drag)
      window.removeEventListener(`touchmove`, handle_drag)
      window.removeEventListener(`mouseup`, end_drag)
      window.removeEventListener(`touchend`, end_drag)
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
    const dist_p12 = Math.pow(p12.x - pt.x, 2) + Math.pow(p12.y - pt.y, 2)
    const dist_p23 = Math.pow(p23.x - pt.x, 2) + Math.pow(p23.y - pt.y, 2)
    const dist_p31 = Math.pow(p31.x - pt.x, 2) + Math.pow(p31.y - pt.y, 2)

    if (dist_p12 <= dist_p23 && dist_p12 <= dist_p31) {
      return p12
    } else if (dist_p23 <= dist_p12 && dist_p23 <= dist_p31) {
      return p23
    } else {
      return p31
    }
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
  <span class="metric-name">
    {METRICS.CPS.label}
    <Tooltip tip_style="z-index: 20; font-size: 0.8em;">
      <svg style="opacity: 0.7; cursor: help;"><use href="#icon-info" /></svg>
      {#snippet tip()}
        {@html METRICS.CPS.description}
      {/snippet}
    </Tooltip>
  </span>

  <button class="reset-button" onclick={reset_weights} title="Reset to default weights">
    Reset
  </button>

  <!-- svelte-ignore a11y_click_events_have_key_events -->
  <!-- svelte-ignore a11y_no_noninteractive_element_interactions -->
  <svg
    bind:this={svg_element}
    width={size}
    height={size}
    viewBox="0 0 {size} {size}"
    onclick={handle_svg_click}
    role="img"
    aria-label="Radar chart for adjusting metric weights"
  >
    <!-- Axes -->
    {#each Object.values(CPS_CONFIG) as weight, idx (weight.label)}
      {@const angle = (2 * Math.PI * idx) / Object.values(CPS_CONFIG).length}
      {@const x = center.x + Math.cos(angle) * radius * 0.8}
      {@const y = center.y + Math.sin(angle) * radius * 0.8}
      {@const label_radius = idx === 2 ? radius * 1.0 : radius * 0.9}
      {@const label_x = center.x + Math.cos(angle) * label_radius}
      {@const label_y = center.y + Math.sin(angle) * label_radius}
      {@const spacing = idx === 1 ? `1.5em` : `1.2em`}

      <line
        x1={center.x}
        y1={center.y}
        x2={x}
        y2={y}
        stroke="rgba(255, 255, 255, 0.4)"
        stroke-width="1"
      />

      <!-- Axis labels -->
      <text
        x={label_x}
        y={label_y}
        text-anchor="middle"
        dominant-baseline="middle"
        font-size="14"
        fill={colors[idx]}
      >
        <!-- Handle subscripts and superscripts manually since <sub> and <sup> are not supported in SVG -->
        {#if weight.label.includes(`<sub>`)}
          {@const parts = weight.label.split(/<sub>|<\/sub>/)}
          {parts[0]}
          <tspan baseline-shift="sub" font-size="10">{parts[1]}</tspan>
        {:else if weight.label.includes(`<sup>`)}
          {@const parts = weight.label.split(/<sup>|<\/sup>/)}
          {parts[0]}
          <tspan baseline-shift="super" font-size="10">{parts[1]}</tspan>
        {:else}
          {@html weight.label}
        {/if}
        <tspan dy={spacing} x={label_x} font-size="12" font-weight="bold"
          >{((weight.weight as number) * 100).toFixed(0)}%</tspan
        >
      </text>
    {/each}

    <!-- Triangle area -->
    {#if Object.values(CPS_CONFIG).length === 3}
      <path
        d="M {axis_points[0].x} {axis_points[0].y} L {axis_points[1].x} {axis_points[1]
          .y} L {axis_points[2].x} {axis_points[2].y} Z"
        fill="rgba(255, 255, 255, 0.1)"
        stroke="rgba(255, 255, 255, 0.3)"
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
        stroke="rgba(255, 255, 255, 0.1)"
        stroke-width="1"
      />
    {/each}

    <!-- Colored areas for each metric -->
    {#if Object.values(CPS_CONFIG).length === 3}
      {@const [{ x, y }, { x: px, y: py }] = [center, point]}
      {#each axis_points as { x: x0, y: y0 }, idx ([x0, y0])}
        {@const path = `M ${x} ${y} L ${x0} ${y0} L ${px} ${py} Z`}
        <path d={path} fill={colors[idx]} stroke="none" opacity="0.5" />
      {/each}
    {/if}

    <!-- Draggable knob: first element is larger invisible hit area for the smaller visible knob above it -->
    {#each [{ fill: `transparent`, r: 20 }, { fill: `white`, stroke: `black`, r: 8 }] as knob_style (knob_style.r)}
      <circle
        cx={point.x}
        cy={point.y}
        style="cursor: move;"
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
    padding: 1em 1em 0 0;
    margin: 0;
    position: relative;
    background: var(--light-bg);
    border-radius: 4px;
  }
  svg {
    touch-action: none; /* Prevents default touch behaviors */
    cursor: pointer; /* Show pointer cursor to hint clickability */
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
    border: 1px solid rgba(255, 255, 255, 0.15);
    border-radius: 3px;
    cursor: pointer;
    padding: 0.15em 0.35em;
    font-size: 0.8em;
    margin-left: auto;
  }
  .reset-button:hover {
    background: rgba(255, 255, 255, 0.05);
  }
</style>
