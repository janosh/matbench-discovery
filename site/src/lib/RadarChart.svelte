<script lang="ts">
  import { onMount } from 'svelte'
  import type { MetricWeight } from './types'

  // Define props interface
  interface Props {
    weights: MetricWeight[]
    size?: number
    onchange?: (new_weights: MetricWeight[]) => void | undefined
  }
  let { weights = [], size = 200, onchange = undefined }: Props = $props()

  // State for the draggable point
  let is_dragging = $state(false)
  let point = $state({ x: 0, y: 0 })
  let svg_element: SVGSVGElement
  let radius = size / 2
  let center = { x: radius, y: radius }

  const colors = [
    `rgba(255, 99, 132, 0.7)`, // red for F1
    `rgba(255, 206, 86, 0.7)`, // yellow for kappa
    `rgba(54, 162, 235, 0.7)`, // blue for RMSD
  ]

  // Compute axes points coordinates
  let axis_points = $derived(
    weights.map((_, idx) => {
      const angle = (2 * Math.PI * idx) / weights.length
      return {
        x: center.x + Math.cos(angle) * radius * 0.8,
        y: center.y + Math.sin(angle) * radius * 0.8,
      }
    }),
  )

  // Initialize point position from weights
  $effect(() => {
    update_point_from_weights(weights)
  })

  function update_point_from_weights(current_weights: MetricWeight[]) {
    if (!current_weights || current_weights.length < 3) return

    // For 3 axes, we can use barycentric coordinates
    const points = axis_points

    // Calculate weighted position
    let { x, y } = center
    for (let idx = 0; idx < current_weights.length; idx++) {
      const weight = current_weights[idx].value
      x += (points[idx].x - center.x) * weight
      y += (points[idx].y - center.y) * weight
    }
    // Update the point with new coordinates
    point = { x, y }
  }

  function update_weights_from_point() {
    // Calculate weights using barycentric coordinates for triangular space
    if (weights.length !== 3) {
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

    // Update weights
    const new_weights = weights.map((wt, idx) => ({
      ...wt,
      value: new_values[idx],
    }))

    // Notify parent component
    onchange?.(new_weights)
  }

  // Helper to calculate triangle area using cross product
  function calc_triangle_area(
    p1: { x: number; y: number },
    p2: { x: number; y: number },
    p3: { x: number; y: number },
  ) {
    return Math.abs(
      (p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y)) / 2,
    )
  }

  // Handle click on SVG to jump to position
  function handle_svg_click(event: MouseEvent) {
    // Get click coordinates relative to SVG
    const target = event.currentTarget as SVGSVGElement
    const svg_rect = target.getBoundingClientRect()
    const click_x = event.clientX - svg_rect.left
    const click_y = event.clientY - svg_rect.top

    // Check if click is inside the triangle
    const [a, b, c] = axis_points
    if (is_point_in_triangle({ x: click_x, y: click_y }, a, b, c)) {
      // Update the draggable point position
      point = { x: click_x, y: click_y }

      // on click, immediately update weights based on new position
      update_weights_from_point()
    }
  }

  // Move point to a position with triangle constraints
  // during dragging, don't update weights - only move the point visually
  // this prevents table rerendering during drag which causes the viewport to scroll (terrible UX)
  function move_point_to_position(x: number, y: number) {
    const [a, b, c] = axis_points
    if (weights.length === 3 && is_point_in_triangle({ x, y }, a, b, c)) {
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

  // Helper to check if a point is inside a triangle
  function is_point_in_triangle(
    p: { x: number; y: number },
    a: { x: number; y: number },
    b: { x: number; y: number },
    c: { x: number; y: number },
  ) {
    // Compute barycentric coordinates
    const denominator = (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y)
    const alpha = ((b.y - c.y) * (p.x - c.x) + (c.x - b.x) * (p.y - c.y)) / denominator
    const beta = ((c.y - a.y) * (p.x - c.x) + (a.x - c.x) * (p.y - c.y)) / denominator
    const gamma = 1 - alpha - beta

    // If all coordinates are between 0 and 1, point is inside
    return alpha >= 0 && beta >= 0 && gamma >= 0 && alpha <= 1 && beta <= 1 && gamma <= 1
  }

  // Helper to find closest point on triangle
  function get_closest_point_on_triangle(
    p: { x: number; y: number },
    a: { x: number; y: number },
    b: { x: number; y: number },
    c: { x: number; y: number },
  ) {
    // First, find the closest point on each edge
    const ab = closest_point_on_line(p, a, b)
    const bc = closest_point_on_line(p, b, c)
    const ca = closest_point_on_line(p, c, a)

    // Then, find which of those points is closest to p
    const dist_ab = Math.pow(ab.x - p.x, 2) + Math.pow(ab.y - p.y, 2)
    const dist_bc = Math.pow(bc.x - p.x, 2) + Math.pow(bc.y - p.y, 2)
    const dist_ca = Math.pow(ca.x - p.x, 2) + Math.pow(ca.y - p.y, 2)

    if (dist_ab <= dist_bc && dist_ab <= dist_ca) {
      return ab
    } else if (dist_bc <= dist_ab && dist_bc <= dist_ca) {
      return bc
    } else {
      return ca
    }
  }

  // Helper to find closest point on a line segment
  function closest_point_on_line(
    p: { x: number; y: number },
    a: { x: number; y: number },
    b: { x: number; y: number },
  ) {
    const atob = { x: b.x - a.x, y: b.y - a.y }
    const atop = { x: p.x - a.x, y: p.y - a.y }
    const len = atob.x * atob.x + atob.y * atob.y
    let dot = atop.x * atob.x + atop.y * atob.y
    const t = Math.max(0, Math.min(1, dot / len))

    return {
      x: a.x + atob.x * t,
      y: a.y + atob.y * t,
    }
  }

  function end_drag() {
    is_dragging = false
    window.removeEventListener(`mousemove`, handle_drag)
    window.removeEventListener(`touchmove`, handle_drag)
    window.removeEventListener(`mouseup`, end_drag)
    window.removeEventListener(`touchend`, end_drag)
  }

  onMount(() => {
    update_point_from_weights(weights)
  })
</script>

<div class="radar-chart-container">
  <!-- svelte-ignore a11y_click_events_have_key_events -->
  <!-- svelte-ignore a11y_no_noninteractive_element_interactions -->
  <svg
    bind:this={svg_element}
    width={size}
    height={size}
    viewBox={`0 0 ${size} ${size}`}
    onclick={handle_svg_click}
    role="img"
    aria-label="Radar chart for adjusting metric weights"
  >
    <!-- Axes -->
    {#each weights as weight, idx}
      {@const angle = (2 * Math.PI * idx) / weights.length}
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
          >{(weight.value * 100).toFixed(0)}%</tspan
        >
      </text>
    {/each}

    <!-- Triangle area -->
    {#if weights.length === 3}
      <path
        d={`M ${axis_points[0].x} ${axis_points[0].y} L ${axis_points[1].x} ${axis_points[1].y} L ${axis_points[2].x} ${axis_points[2].y} Z`}
        fill="rgba(255, 255, 255, 0.1)"
        stroke="rgba(255, 255, 255, 0.3)"
        stroke-width="1"
      />
    {/if}

    <!-- Background circular grid -->
    {#each [0.2, 0.4, 0.6, 0.8] as grid_radius}
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
    {#if weights.length === 3}
      {@const { x, y } = center}
      {@const { x: x0, y: y0 } = axis_points[0]}
      {@const { x: x1, y: y1 } = axis_points[1]}
      {@const { x: x2, y: y2 } = axis_points[2]}
      <path
        d={`M ${x} ${y} L ${x0} ${y0} L ${point.x} ${point.y} Z`}
        fill={colors[0]}
        stroke="none"
        opacity="0.5"
      />
      <path
        d="M {x} {y} L {x1} {y1} L {point.x} {point.y} Z"
        fill={colors[1]}
        stroke="none"
        opacity="0.5"
      />
      <path
        d="M {x} {y} L {x2} {y2} L {point.x} {point.y} Z"
        fill={colors[2]}
        stroke="none"
        opacity="0.5"
      />
    {/if}

    <!-- svelte-ignore a11y_interactive_supports_focus -->
    <!-- Draggable knob -->
    <circle
      cx={point.x}
      cy={point.y}
      r="8"
      fill="white"
      stroke="black"
      stroke-width="2"
      style="cursor: move;"
      onmousedown={start_drag}
      ontouchstart={start_drag}
      role="button"
      aria-label="Drag to adjust weight balance"
    />
  </svg>
</div>

<style>
  .radar-chart-container {
    display: flex;
    flex-direction: column;
    align-items: center;
    padding: 0.2em;
    margin: 0;
  }
  svg {
    touch-action: none; /* Prevents default touch behaviors */
    cursor: pointer; /* Show pointer cursor to hint clickability */
  }
</style>
