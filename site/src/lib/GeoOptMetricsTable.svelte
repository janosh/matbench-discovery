<script lang="ts">
  import { HeatmapTable } from '$lib'
  import type { ModelData } from '$lib/types'
  import { ALL_METRICS, GEO_OPT_SYMMETRY_METRICS, METADATA_COLS } from './labels'
  import { assemble_row_data } from './metrics'

  let { column_order = $bindable([]), ...rest } = $props()
  let columns = $derived([
    METADATA_COLS.model_name,
    ALL_METRICS.RMSD,
    ...Object.values(GEO_OPT_SYMMETRY_METRICS).map((col) => ({
      ...col,
      visible: true,
    })),
  ])

  const discovery_set = `full_test_set`
  const model_filter = (_model: ModelData) => true
  const show_energy_only = $state(false)
  const show_non_compliant = $state(true)
  const show_compliant = $state(true)

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
  {...rest}
  style="margin: 2em 0"
/>
