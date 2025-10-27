<script lang="ts">
  import { HeatmapTable, TableControls } from '$lib'
  import type { Label, ModelData } from '$lib/types'
  import type { HTMLAttributes } from 'svelte/elements'
  import {
    ALL_METRICS,
    GEO_OPT_SYMMETRY_METRICS,
    HYPERPARAMS,
    METADATA_COLS,
  } from './labels'
  import { assemble_row_data } from './metrics'

  let {
    column_order = $bindable([]),
    show_heatmap = $bindable(true),
    show_compliant = $bindable(true),
    show_non_compliant = $bindable(true),
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    column_order?: string[]
    show_heatmap?: boolean
    show_compliant?: boolean
    show_non_compliant?: boolean
  } = $props()

  const hyperparam_cols = [
    HYPERPARAMS.ase_optimizer,
    HYPERPARAMS.max_steps,
    HYPERPARAMS.max_force,
    HYPERPARAMS.cell_filter,
  ]
  let columns = $state<Label[]>([
    METADATA_COLS.model_name,
    ALL_METRICS.RMSD,
    ...Object.values(GEO_OPT_SYMMETRY_METRICS).map((col) => ({
      ...col,
      visible: true,
    })),
    ...hyperparam_cols.map((col) => ({ ...col, visible: true, sortable: true })),
  ])

  const discovery_set = `full_test_set`
  const model_filter = (_model: ModelData) => true
  const show_energy_only = false

  // recalculate metrics_data whenever filter settings change
  let metrics_data = $derived(
    assemble_row_data(
      discovery_set,
      model_filter,
      show_energy_only,
      show_non_compliant,
      show_compliant,
    ),
  )
</script>

<HeatmapTable
  data={metrics_data}
  {columns}
  bind:column_order
  bind:show_heatmap
  {...rest}
>
  {#snippet controls()}
    <TableControls
      bind:columns
      bind:show_heatmap
      bind:show_compliant
      bind:show_non_compliant
    />
  {/snippet}
</HeatmapTable>
