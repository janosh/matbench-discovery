<script lang="ts">
  import MetricsTable from '$lib/table/MetricsTable.svelte'
  import type { Label } from '$lib/types'
  import type { ComponentProps } from 'svelte'
  import {
    ALL_METRICS,
    GEO_OPT_SYMMETRY_METRICS,
    HYPERPARAMS,
    METADATA_COLS,
  } from '../labels'

  let {
    column_order = $bindable([]),
    sort = $bindable({ column: ALL_METRICS.RMSD.key, dir: `asc` }),
    ...rest
  }: ComponentProps<typeof MetricsTable> = $props()

  const columns: Label[] = [
    METADATA_COLS.model_name,
    ALL_METRICS.RMSD,
    ...Object.values(GEO_OPT_SYMMETRY_METRICS).map((col) => ({
      ...col,
      group: `Symmetry`,
      visible: true,
    })),
    ...[
      HYPERPARAMS.ase_optimizer,
      HYPERPARAMS.max_steps,
      HYPERPARAMS.max_force,
      HYPERPARAMS.cell_filter,
      HYPERPARAMS.n_layers,
      HYPERPARAMS.graph_construction_radius,
    ].map((col, idx) => ({
      ...col,
      group: `Hyperparams`,
      sortable: true,
      visible: idx < 4,
    })),
  ].map((col) => ({
    ...col,
    // Geometry hyperparameters show their units in the header.
    label: col.unit
      ? `${col.label} <span style="font-weight: 200">(${col.unit})</span>`
      : col.label,
  }))
</script>

<MetricsTable
  {...rest}
  discovery_set="full_test_set"
  model_filter={(model) => model.metrics?.geo_opt != null}
  column_labels={columns}
  show_row_numbers={false}
  bind:sort
  bind:column_order
/>
