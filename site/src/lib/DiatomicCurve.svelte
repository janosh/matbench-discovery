<script lang="ts">
  import { MODELS } from '$lib'
  import { ScatterPlot, type InternalPoint } from 'matterviz'

  interface Props {
    formula: string
    curves: Array<{
      model_key: string
      distances: number[]
      energies: number[]
      color: string
    }>
    tooltip_point?: InternalPoint | null
    hovered?: boolean
    [key: string]: unknown
  }
  let {
    formula,
    curves,
    tooltip_point = $bindable(null),
    hovered = $bindable(false),
    ...rest
  }: Props = $props()

  // Create a map of model keys to labels from MODELS
  // Extract base model name by removing version suffix
  const model_labels = new Map(
    MODELS.map((model) => {
      const base_name = model.dirname
      const label = model.model_name ?? base_name
      // Map both the base name and any known version-specific names to the same label
      return [base_name, label]
    }),
  )

  // Function to get model label, handling version suffixes
  function get_model_label(model_key: string): string {
    // Try to match the model key to known base names
    const base_name = MODELS.find((model) =>
      model_key?.startsWith(model.dirname),
    )?.dirname

    return base_name ? (model_labels.get(base_name) ?? model_key) : model_key
  }

  const x_lim: [number, number] = [0.2, 6]
  const y_lim: [number, number] = [-8, 20]

  let series = $derived(
    curves.map((curve) => {
      // Filter out points outside the x_lim range
      const filtered_indices = curve.distances
        .map((dist, idx) => [dist, idx])
        .filter(([dist]) => dist >= x_lim[0] && dist <= x_lim[1])
        .map(([_, idx]) => idx)

      const filtered_distances = filtered_indices.map((i) => curve.distances[i])
      const filtered_energies = filtered_indices.map((i) => curve.energies[i])

      // Shift energies so the energy at infinite separation (last point) is 0
      const shifted_energies = filtered_energies.map(
        (energy) => energy - filtered_energies[filtered_energies.length - 1],
      )

      return {
        x: filtered_distances,
        y: shifted_energies,
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

<!-- TODO increase font size of axes titles and tick labels -->
<div class="plot" {...rest}>
  <h3>{formula}</h3>
  <ScatterPlot
    {series}
    x_label="Distance (Å)"
    y_label="Energy (eV)"
    x_format=".1f"
    y_format=".2f"
    markers="line+points"
    {x_lim}
    {y_lim}
    bind:tooltip_point
    bind:hovered
    legend={null}
  >
    {#snippet tooltip({ x_formatted, y_formatted, metadata })}
      <div
        style="min-width: 10em; background: rgba(255, 255, 255, 0.1); padding: 2pt 4pt; border-radius: 3pt;"
      >
        <strong>{metadata?.model_label ?? ``}</strong>
        <br />
        Distance = {x_formatted} Å
        <br />
        Energy = {y_formatted} eV
      </div>
    {/snippet}
  </ScatterPlot>
</div>

<style>
  .plot {
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
