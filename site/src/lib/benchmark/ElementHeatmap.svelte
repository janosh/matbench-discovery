<script lang="ts">
  import type { BenchmarkDataset } from './data'
  import { format_num } from 'matterviz/labels'
  import { PeriodicTable } from 'matterviz/periodic-table'

  let { dataset }: { dataset: BenchmarkDataset } = $props()
  const counts = $derived([
    ...new Set(Object.values(dataset.element_counts).filter((count) => count !== null)),
  ])
  const log_scale = $derived(dataset.id !== `diatomics` && counts.length > 1)
</script>

<PeriodicTable
  heatmap_values={dataset.element_counts}
  log={log_scale}
  lanth_act_tiles={[]}
  tile_props={{ show_name: false, show_number: false }}
  missing={{ color: `var(--text-color)`, style: `opacity: 0.15` }}
  color_bar_props={{
    title: `${dataset.count_unit}${log_scale ? ` (log)` : ``}`,
    tick_labels: counts.length === 1 ? counts : 3,
    snap_ticks: false,
    tick_format: `d`,
  }}
  style="--ptable-min-tile-size: 0; --elem-symbol-font-size: 60cqw"
>
  {#snippet tooltip({ element, value })}
    <strong>{element.name}</strong><br />
    {dataset.count_unit} containing {element.symbol}: {format_num(
      Number(value ?? 0),
      `,`,
    )}
  {/snippet}
</PeriodicTable>
