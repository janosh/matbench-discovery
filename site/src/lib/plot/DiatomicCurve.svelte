<script lang="ts">
  import { MODELS } from '$lib'
  import { type InternalPoint, ScatterPlot } from 'matterviz'
  import { element_data, type ChemicalElement } from 'matterviz/element'
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
      line_width?: number
    }[]
    tooltip_point?: InternalPoint | null
    hovered?: boolean
  } = $props()

  function get_model_label(model_key: string): string {
    const model =
      MODELS.find(
        (model_entry) =>
          model_entry.model_name === model_key || model_entry.dirname === model_key,
      ) ??
      MODELS.toSorted(
        (model_a, model_b) => model_b.dirname.length - model_a.dirname.length,
      ).find((model_entry) => model_key.startsWith(model_entry.dirname))
    return model ? (model.model_name ?? model.dirname) : model_key
  }

  const x_range: [number, number] = [0.2, 6]
  const y_range: [number, number] = [-8, 20]
  const no_tween = { duration: 0 } as const
  const elements = element_data as ChemicalElement[]
  let element_symbol = $derived(formula.split(`-`, 1)[0])
  let element = $derived(
    elements.find((element_entry) => element_entry.symbol === element_symbol),
  )
  let header_text = $derived(element ? `${element.number} ${formula}` : formula)
  let header_title = $derived(
    element
      ? `${element.name} (Z=${element.number})\n` +
          `Row=${element.row} Col=${element.column}`
      : formula,
  )

  let series = $derived(
    curves.map((curve) => {
      const model_label = get_model_label(curve.model_key)
      // Keep only points within the x range, pairing each distance with its energy
      const points = curve.distances
        .map((distance, idx) => ({ distance, energy: curve.energies[idx] }))
        .filter(({ distance }) => distance >= x_range[0] && distance <= x_range[1])
      // Shift energies so the energy at infinite separation (last point) is 0
      const ref_energy = points.at(-1)?.energy ?? 0

      return {
        id: curve.model_key,
        label: model_label,
        x: points.map((point) => point.distance),
        y: points.map((point) => point.energy - ref_energy),
        markers: `line+points` as const,
        metadata: {
          model_key: curve.model_key,
          model_label,
        },
        point_style: { fill: curve.color, radius: 1.5, stroke_width: 0 },
        line_style: {
          stroke: curve.color,
          // DFT references (line_width set) get a thicker line to stand out
          ...(curve.line_width ? { stroke_width: curve.line_width } : {}),
        },
      }
    }),
  )
</script>

<div {...rest} class="plot {rest.class ?? ``}">
  <h3 aria-label={header_title}>
    {header_text}
    <span class="element-tooltip" role="tooltip">{header_title}</span>
  </h3>
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
    position: relative;
    align-self: center;
    margin: 0 0 -5pt;
    text-align: center;
    font-size: 0.9em;
  }
  .element-tooltip {
    position: absolute;
    z-index: 10;
    bottom: calc(100% + 0.4em);
    left: 50%;
    padding: 0.35em 0.55em;
    border: 1px solid
      var(--border-color, color-mix(in srgb, currentColor 25%, transparent));
    border-radius: 4px;
    background: var(--page-bg, var(--card-bg, white));
    box-shadow: 0 2px 8px rgb(0 0 0 / 15%);
    color: var(--text-color, currentColor);
    font-size: 0.85rem;
    font-weight: 400;
    line-height: 1.3;
    opacity: 0;
    pointer-events: none;
    transform: translateX(-50%);
    transition: opacity 0.12s ease;
    visibility: hidden;
    white-space: pre;
  }
  h3:hover .element-tooltip {
    opacity: 1;
    visibility: visible;
  }
</style>
