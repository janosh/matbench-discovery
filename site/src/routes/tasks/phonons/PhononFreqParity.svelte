<script lang="ts">
  // Quantile parity between ML and DFT phonon frequency spectra for one model.
  import type kappa_data from '$figs/kappa-103-analysis.jsonl'
  import { parity_diagonal } from '$lib/fig-helpers'
  import { format_num, sanitize_compact_formula } from 'matterviz'
  import { DEFAULT_SERIES_SYMBOLS, ScatterPlot } from 'matterviz/plot'
  import type { DataSeries } from 'matterviz/plot'
  import { spacegroup_num_to_crystal_sys } from 'matterviz/symmetry'
  import type { HTMLAttributes } from 'svelte/elements'

  type KappaModelEntry = (typeof kappa_data)[`models`][number]
  const n_quantiles = 17

  interface FreqParityPoint extends Record<string, unknown> {
    material_id: string
    formula: string
    spg_num: number
    quantile: number
    freq_w1: number
    crystal_system: string
  }
  let {
    entry,
    base,
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    entry: KappaModelEntry
    base: Pick<typeof kappa_data, `material_ids` | `formulas` | `spg_nums`>
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
  let points = $derived.by((): FreqParityPoint[] => {
    const result: FreqParityPoint[] = []
    let pair_idx = 0
    for (const [mat_idx, material_id] of base.material_ids.entries()) {
      const freq_w1 = entry.freq_w1[mat_idx]
      if (freq_w1 === null) continue
      const spg_num = base.spg_nums[mat_idx]
      const system = spacegroup_num_to_crystal_sys(spg_num)
      const dft = pairs.dft.slice(pair_idx, pair_idx + n_quantiles)
      const ml = pairs.ml.slice(pair_idx, pair_idx + n_quantiles)
      pair_idx += n_quantiles
      if (!system || dft.length !== n_quantiles || ml.length !== n_quantiles) continue

      result.push(
        ...dft.map((_, quantile_idx) => ({
          material_id,
          formula: base.formulas[mat_idx],
          spg_num,
          quantile: quantile_idx / (n_quantiles - 1),
          freq_w1,
          crystal_system: system,
        })),
      )
    }
    return result
  })
  let series = $derived<DataSeries<FreqParityPoint>[]>([
    {
      x: pairs.dft,
      y: pairs.ml,
      metadata: points,
      label: entry.label,
      markers: `points` as const,
      color_values: points.map((point) => point.freq_w1),
      point_style: points.map((point) => ({
        radius: 3,
        stroke: `white`,
        stroke_width: 0.4,
        symbol_type:
          DEFAULT_SERIES_SYMBOLS[
            Math.floor(point.quantile * 3) % DEFAULT_SERIES_SYMBOLS.length
          ],
      })),
    },
  ])
</script>

{#if pairs.dft.length}
  <ScatterPlot
    {series}
    ref_lines={[parity_diagonal]}
    x_axis={{ label: `PBE phonon freq. (THz)`, range: extent, format: `.3~` }}
    y_axis={{ label: `${entry.label} phonon freq. (THz)`, range: extent, format: `.3~` }}
    color_bar={{ title: `W₁(ω) (THz)`, tick_format: `.3~` }}
    {...rest}
  >
    {#snippet tooltip({ x_formatted, y_formatted, metadata })}
      {#if metadata}
        {@const point = metadata as FreqParityPoint}
        <strong>{point.material_id}</strong>
        {@html sanitize_compact_formula(point.formula)} ({point.crystal_system}, SG
        {point.spg_num})<br />
        quantile: {format_num(point.quantile, `.0%`)}<br />
        spectrum W1: {format_num(point.freq_w1, `.3~`)} <small>THz</small><br />
      {/if}
      PBE: {x_formatted} THz<br />
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
