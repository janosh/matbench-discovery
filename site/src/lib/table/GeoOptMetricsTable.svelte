<script lang="ts">
  import { CompareToggle, ModelRowMenu, TableControls } from '$lib'
  import { mark_compared_rows, toggle_row_model } from '$lib/model-comparison.svelte'
  import { make_table_filters } from '$lib/models.svelte'
  import { METRICS_TABLE_ROOT_STYLE } from '$lib/table/MetricsTable.svelte'
  import type { SortState, UrlTableFilters } from '$lib/url-state.svelte'
  import type { Label, TableLabel } from '$lib/types'
  import { HeatmapTable, type CellSnippetArgs, type RowData } from 'matterviz'
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
    sort = $bindable({ column: ALL_METRICS.RMSD.key, dir: `asc` }),
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    column_order?: string[]
    filters?: UrlTableFilters
    sort?: SortState
  } = $props()
  // toggled from TableControls; no page binds it, so plain local state
  let show_selected_only = $state(false)

  // Append unit in thin font and (higher/lower=better) hint to column tooltip
  function enrich_col(col: Label, overrides: Partial<Label> = {}): TableLabel {
    let { label } = col
    if (col.unit) {
      label = `${label} <span style="font-weight: 200">(${col.unit})</span>`
    }
    const better = overrides.better ?? col.better ?? undefined
    const description = append_better_hint(col)
    return { ...col, ...overrides, better, description, label }
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

  let columns = $state<TableLabel[]>([
    enrich_col(METADATA_COLS.model_name),
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

  let metrics_data = $derived(
    mark_compared_rows(
      assemble_row_data(
        `full_test_set`,
        (model) => model.metrics?.geo_opt != null,
        filters.matches,
      ),
      show_selected_only,
    ).map((row) => {
      for (const [from, to] of Object.entries(key_remap)) {
        if (from in row) row[to] = row[from]
      }
      return row
    }),
  )
</script>

<!-- model link plus a compare button that appears on row hover -->
{#snippet model_cell({ row, val }: CellSnippetArgs)}
  {@const { model_key, model_name } = row as { model_key: string; model_name: string }}
  {@html String(val)}
  <CompareToggle {model_key} {model_name} compact />
{/snippet}

<ModelRowMenu>
  <HeatmapTable
    data={metrics_data as RowData[]}
    {columns}
    bind:sort
    special_cells={{ [METADATA_COLS.model_name.label]: model_cell }}
    default_num_format=".3f"
    bind:column_order
    bind:show_heatmap={filters.show_heatmap}
    on_row_double_click={toggle_row_model}
    {...rest}
    class={[`leaderboard`, rest.class]}
    root_style={METRICS_TABLE_ROOT_STYLE}
  >
    {#snippet controls()}
      <!-- z-index > 2 to sit above sticky table headers (z-index: 2) -->
      <TableControls
        bind:columns
        bind:show_selected_only
        {filters}
        style="position: relative; z-index: 5"
      />
    {/snippet}
  </HeatmapTable>
</ModelRowMenu>
