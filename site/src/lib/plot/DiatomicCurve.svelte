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

  function get_model_label(model_key: string): string {
    const model = MODELS.find((mdl) => model_key.startsWith(mdl.dirname))
    return model ? (model.model_name ?? model.dirname) : model_key
  }

  const x_range: [number, number] = [0.2, 6]
  const y_range: [number, number] = [-8, 20]
  const no_tween = { duration: 0 } as const

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
        point_style: { radius: 1.5, stroke_width: 0 },
      }
    }),
  )
</script>

<div {...rest} class="plot {rest.class ?? ``}">
  <h3>{formula}</h3>
  <ScatterPlot
    {series}
    x_axis={{
      label: `Distance (Å)`,
      format: `.1f`,
      range: x_range,
      label_shift: { y: -30 },
    }}
    y_axis={{ label: `Energy (eV)`, format: `.2f`, range: y_range }}
    bind:tooltip_point
    bind:hovered
    legend={null}
    point_tween={no_tween}
    line_tween={no_tween}
  >
    {#snippet tooltip({ x_formatted, y_formatted, metadata })}
      <strong>{metadata?.model_label ?? ``}</strong><br />
      Distance = {x_formatted} Å<br />
      Energy = {y_formatted} eV
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
  }
  h3 {
    margin: 0 0 -5pt;
    text-align: center;
    font-size: 0.9em;
  }
</style>
