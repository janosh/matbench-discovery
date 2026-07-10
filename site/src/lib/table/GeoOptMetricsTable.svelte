<script lang="ts">
  import { TableControls } from '$lib'
  import { make_table_filters } from '$lib/models.svelte'
  import { TABLE_STYLE_VARS } from '$lib/table/MetricsTable.svelte'
  import type { UrlTableFilters } from '$lib/url-state.svelte'
  import type { Label, ModelData } from '$lib/types'
  import { HeatmapTable, type Label as MattervizLabel } from 'matterviz'
  import type { HTMLAttributes } from 'svelte/elements'
  import {
    ALL_METRICS,
    GEO_OPT_SYMMETRY_METRICS,
    HYPERPARAMS,
    METADATA_COLS,
  } from '../labels'
  import { append_better_hint, assemble_row_data } from '../metrics'

  let {
    column_order = $bindable([]),
    filters = make_table_filters(),
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    column_order?: string[]
    filters?: UrlTableFilters
  } = $props()

  // Append unit in thin font and (higher/lower=better) hint to column tooltip
  function enrich_col(col: Label, overrides: Partial<Label> = {}): Label {
    let { label } = col
    if (col.unit) {
      label = `${label} <span style="font-weight: 200">(${col.unit})</span>`
    }
    return { ...col, ...overrides, label, description: append_better_hint(col) }
  }

  // Define grouped columns: [source_cols, group_name, extra_overrides]
  const grouped_defs: [Label[], string, Partial<Label>][] = [
    [Object.values(GEO_OPT_SYMMETRY_METRICS), `Symmetry`, { visible: true }],
    [
      [
        HYPERPARAMS.ase_optimizer,
        HYPERPARAMS.max_steps,
        HYPERPARAMS.max_force,
        HYPERPARAMS.cell_filter,
        HYPERPARAMS.n_layers,
        HYPERPARAMS.graph_construction_radius,
      ],
      `Hyperparams`,
      { sortable: true },
    ],
  ]

  // First 4 hyperparams visible by default, last 2 hidden
  const hidden_hyperparam_keys = new Set([
    HYPERPARAMS.n_layers.key,
    HYPERPARAMS.graph_construction_radius.key,
  ])

  let columns = $state<Label[]>([
    METADATA_COLS.model_name,
    enrich_col(ALL_METRICS.RMSD),
    ...grouped_defs.flatMap(([cols, group, extras]) =>
      cols.map((col) =>
        enrich_col(col, {
          ...extras,
          visible: extras.visible ?? !hidden_hyperparam_keys.has(col.key),
          group,
        }),
      ),
    ),
  ])

  // HeatmapTable uses `"${key} (${group})"` as col ID when group is set,
  // So remap data keys to match
  const key_remap: Record<string, string> = Object.fromEntries(
    grouped_defs.flatMap(([cols, group]) =>
      cols.map((col) => [col.key, `${col.key} (${group})`]),
    ),
  )

  const has_geo_opt_metrics = (model: ModelData): boolean =>
    model.metrics?.geo_opt != null && typeof model.metrics.geo_opt === `object`

  let metrics_data = $derived(
    assemble_row_data(`full_test_set`, has_geo_opt_metrics, filters.matches).map(
      (row) => {
        for (const [from, to] of Object.entries(key_remap)) {
          if (from in row) row[to] = row[from]
        }
        return row
      },
    ),
  )
</script>

<HeatmapTable
  data={metrics_data}
  columns={columns as MattervizLabel[]}
  initial_sort={{ column: ALL_METRICS.RMSD.key, direction: `asc` }}
  default_num_format=".3f"
  bind:column_order
  bind:show_heatmap={filters.show_heatmap}
  {...rest}
  style="{TABLE_STYLE_VARS}{rest.style ?? ``}"
>
  {#snippet controls()}
    <!-- z-index > 2 to sit above sticky table headers (z-index: 2) -->
    <TableControls bind:columns {filters} style="position: relative; z-index: 5" />
  {/snippet}
</HeatmapTable>
