<script lang="ts">
  // Quantile parity between ML and DFT phonon frequency spectra for one model.
  import type kappa_data from '$figs/kappa-103-analysis.json.gz'
  import { parity_diagonal } from '$lib/fig-helpers'
  import { format_num } from 'matterviz'
  import { ScatterPlot } from 'matterviz/plot'
  import type { DataSeries } from 'matterviz/plot'
  import type { HTMLAttributes } from 'svelte/elements'

  type KappaModelEntry = (typeof kappa_data)[`models`][number]

  let { entry, ...rest }: HTMLAttributes<HTMLDivElement> & {
    entry: KappaModelEntry
  } = $props()

  let pairs = $derived(entry.freq_pairs)
  // Shared range keeps the y=x diagonal corner-to-corner.
  let extent = $derived.by((): [number, number] => {
    if (!pairs.dft.length) return [0, 1]
    const values = [...pairs.dft, ...pairs.ml]
    const [min, max] = [Math.min(...values), Math.max(...values)]
    const pad = (max - min) * 0.05
    return [Math.min(min - pad, 0), max + pad]
  })
  let series = $derived<DataSeries[]>([
    {
      x: pairs.dft,
      y: pairs.ml,
      label: entry.label,
      markers: `points`,
      point_style: { fill: `#4dabf7`, radius: 2.5, stroke: `transparent` },
    },
  ])
</script>

{#if pairs.dft.length}
  <ScatterPlot
    {series}
    ref_lines={[parity_diagonal]}
    x_axis={{ label: `PBE phonon freq. (THz)`, range: extent, format: `.3~` }}
    y_axis={{ label: `${entry.label} phonon freq. (THz)`, range: extent, format: `.3~` }}
    {...rest}
  >
    {#snippet tooltip({ x_formatted, y_formatted })}
      PBE: {x_formatted} THz<br>
      {entry.label}: {y_formatted} THz
    {/snippet}

    {#snippet user_content({ width, height })}
      {#if entry.freq_w1_mean !== null}
        <foreignObject
          x="0"
          y="0"
          {width}
          {height}
          style="pointer-events: none; overflow: visible"
        >
          <div class="plot-annotation">
            spectrum W1 = {format_num(entry.freq_w1_mean, `.3~`)} <small>THz</small>
          </div>
        </foreignObject>
      {/if}
    {/snippet}
  </ScatterPlot>
{:else}
  <p class="plot-state">
    {entry.label} has no phonon frequency data to compare against DFT.
  </p>
{/if}
