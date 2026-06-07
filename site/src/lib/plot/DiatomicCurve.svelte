<script lang="ts">
  import { MODELS } from '$lib'
  import { type InternalPoint, ScatterPlot } from 'matterviz'
  import type { HTMLAttributes } from 'svelte/elements'

  let {
    formula,
    curves,
    tooltip_point = $bindable(null),
    hovered = $bindable(false),
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    formula: string
    curves: {
      model_key: string
      distances: number[]
      energies: number[]
      color: string
    }[]
    tooltip_point?: InternalPoint | null
    hovered?: boolean
  } = $props()

  // Map each model's dirname (base name) to its label; version-specific keys are
  // resolved to a base name in get_model_label below
  const model_labels = new Map(
    MODELS.map((model) => [model.dirname, model.model_name ?? model.dirname]),
  )

  // Function to get model label, handling version suffixes
  function get_model_label(model_key: string): string {
    // Try to match the model key to known base names
    const base_name = MODELS.find((model) => model_key?.startsWith(model.dirname))
      ?.dirname

    return base_name ? (model_labels.get(base_name) ?? model_key) : model_key
  }

  const x_range: [number, number] = [0.2, 6]
  const y_range: [number, number] = [-8, 20]

  let series = $derived(
    curves.map((curve) => {
      // Keep only points within the x range, pairing each distance with its energy
      const points = curve.distances
        .map((distance, idx) => ({ distance, energy: curve.energies[idx] }))
        .filter(({ distance }) => distance >= x_range[0] && distance <= x_range[1])
      // Shift energies so the energy at infinite separation (last point) is 0
      const ref_energy = points.at(-1)?.energy ?? 0

      return {
        x: points.map((point) => point.distance),
        y: points.map((point) => point.energy - ref_energy),
        markers: `line+points` as const,
        metadata: {
          model_key: curve.model_key,
          model_label: get_model_label(curve.model_key),
        },
        point_style: {
          radius: 1.5,
          fill: curve.color,
          fill_opacity: 0.8,
        },
        point_hover: {
          enabled: true,
          scale: 2,
          stroke: `white`,
          stroke_width: 1,
        },
      }
    }),
  )
</script>

<div {...rest} class="plot {rest.class ?? ``}">
  <h3>{formula}</h3>
  <ScatterPlot
    {series}
    x_axis={{ label: `Distance (Å)`, format: `.1f`, range: x_range }}
    y_axis={{ label: `Energy (eV)`, format: `.2f`, range: y_range }}
    bind:tooltip_point
    bind:hovered
    legend={null}
  >
    {#snippet tooltip({ x_formatted, y_formatted, metadata })}
      <div
        style="min-width: 10em; background: rgba(255, 255, 255, 0.1); padding: 2pt 4pt; border-radius: 3pt"
      >
        <strong>{metadata?.model_label ?? ``}</strong><br />
        Distance = {x_formatted} Å<br />
        Energy = {y_formatted} eV
      </div>
    {/snippet}
  </ScatterPlot>
</div>

<style>
  .plot {
    /* bump axis titles + tick labels (matterviz default inherits smaller SVG size) */
    --scatter-font-size: 14px;
    display: flex;
    flex-direction: column;
    box-sizing: border-box;
    height: 300px; /* plot height */
  }
  h3 {
    margin: 0;
    text-align: center;
    font-size: 0.85em;
    font-weight: normal;
  }
</style>
