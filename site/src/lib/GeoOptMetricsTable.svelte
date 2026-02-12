<script lang="ts">
  import { TableControls } from '$lib'
  import type { Label } from '$lib/types'
  import { HeatmapTable, type Label as MattervizLabel } from 'matterviz'
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

  // Append unit in thin font and (higher/lower=better) hint to column tooltip
  function enrich_col(col: Label, overrides: Partial<Label> = {}): Label {
    let { label, description = `` } = col
    if (col.unit) {
      label = `${label} <span style="font-weight: 200">(${col.unit})</span>`
    }
    if (col.better === `higher` || col.better === `lower`) {
      description = description
        ? `${description} (${col.better}=better)`
        : `${col.better}=better`
    }
    return { ...col, label, description, ...overrides }
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
          visible: !hidden_hyperparam_keys.has(col.key),
          ...extras,
          group,
        })
      )
    ),
  ])

  // HeatmapTable uses `"${key} (${group})"` as col ID when group is set,
  // so remap data keys to match
  const key_remap: Record<string, string> = Object.fromEntries(
    grouped_defs.flatMap(([cols, group]) =>
      cols.map((col) => [col.key ?? col.label, `${col.key ?? col.label} (${group})`])
    ),
  )

  // recalculate metrics_data whenever filter settings change
  let metrics_data = $derived(
    assemble_row_data(
      `full_test_set`,
      () => true,
      false,
      show_non_compliant,
      show_compliant,
    )
      .map((row) => {
        const remapped = { ...row }
        for (const [from, to] of Object.entries(key_remap)) {
          if (from in remapped) remapped[to] = remapped[from]
        }
        return remapped
      }),
  )
</script>

<HeatmapTable
  data={metrics_data}
  columns={columns as MattervizLabel[]}
  initial_sort={{ column: `rmsd`, direction: `asc` }}
  sort_hint="Click column headers to sort"
  default_num_format=".3f"
  bind:column_order
  bind:show_heatmap
  {...rest}
>
  {#snippet controls()}
    <!-- z-index > 2 to sit above sticky table headers (z-index: 2) -->
    <div style="position: relative; z-index: 5">
      <TableControls
        bind:columns
        bind:show_heatmap
        bind:show_compliant
        bind:show_non_compliant
      />
    </div>
  {/snippet}
</HeatmapTable>
