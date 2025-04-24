<script lang="ts">
  import { HeatmapTable } from '$lib'
  import type { CellSnippetArgs, Metric, ModelData, RowData } from '$lib/types'
  import { GEO_OPT_SYMMETRY_METRICS, METADATA_COLS, METRICS } from './labels'
  import { assemble_row_data } from './metrics'

  let { ...rest } = $props()

  const links_col: Metric = {
    key: `Links`,
    label: `Links`,
    description: `Links to model resources`,
    style: `border-left: 1px solid black`,
    format: undefined,
    better: null,
    sortable: false,
  }
  // Create HeatmapColumn definitions dynamically
  let columns = $derived([
    METADATA_COLS.model_name,
    METRICS.RMSD,
    ...Object.values(GEO_OPT_SYMMETRY_METRICS),
    links_col,
  ])

  let metrics_data = $state<RowData[]>([])
  const discovery_set = `full_test_set`
  const model_filter = (_model: ModelData) => true
  const show_energy_only = false
  const show_noncompliant = false
  const compliant_clr = `#4caf50`
  const noncompliant_clr = `#4682b4`
  // recalculate metrics_data whenever filter settings, props, or metric_config change
  $effect(() => {
    metrics_data = assemble_row_data(
      discovery_set,
      model_filter,
      show_energy_only,
      show_noncompliant,
      compliant_clr,
      noncompliant_clr,
    )
  })
</script>

{#snippet links_cell({ val }: CellSnippetArgs)}
  {@const links = val as
    | { files: { url: string; title: string; icon: string }[] }
    | null
    | undefined}
  {#if links?.files}
    {#each links.files as { url: href, title, icon }, idx (title + href)}
      {#if href}
        {#if idx > 0}&thinsp;{/if}
        <a {href} {title} target="_blank" rel="noopener noreferrer">
          {@html icon}
        </a>
      {/if}
    {/each}
  {/if}
{/snippet}

<HeatmapTable
  data={metrics_data}
  {columns}
  special_cells={{ Links: links_cell }}
  {...rest}
/>
