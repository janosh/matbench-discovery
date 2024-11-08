<script lang="ts">
  import { HeatmapTable, MODEL_METADATA, model_is_compliant } from '$lib'
  import { geo_opt } from '$root/scripts/metrics-which-is-better.yml'

  export let show_non_compliant: boolean = false
  export let show_metadata: boolean = true
  export let metadata_cols: string[] = []
  export let hide_cols: string[] = []
  export let show_cols = [
    `Model`,
    `RMSD`,
    `Match`,
    `Decrease`,
    `Increase`,
    // `n_structs`,
    ...(show_metadata ? metadata_cols : []),
  ].filter((col) => !hide_cols.includes(col))

  const long_date = (date: string): string =>
    new Date(date).toLocaleDateString(undefined, {
      weekday: `long`,
      year: `numeric`,
      month: `long`,
      day: `numeric`,
    })

  // Transform MODEL_METADATA into table data format
  $: metrics_data = MODEL_METADATA.filter(
    (model) =>
      show_non_compliant || (model_is_compliant(model) && model.metrics?.geo_opt?.rmsd),
  )
    .map((model) => {
      const geo_opt = model.metrics?.geo_opt
      if (!geo_opt) return null

      return {
        Model: `<a title="Version: ${model.model_version}" href="/models/${model.model_key}">${model.model_name}</a>`,
        RMSD: geo_opt.rmsd,
        Match: geo_opt.symmetry_match,
        Decrease: geo_opt.symmetry_decrease,
        Increase: geo_opt.symmetry_increase,
        // n_structs: geo_opt.n_structs,
        'Date Added': `<span title="${long_date(model.date_added)}">${model.date_added}</span>`,
      }
    })
    .sort((row1, row2) => (row1?.RMSD ?? 0) - (row2?.RMSD ?? 0)) // Sort by RMSD ascending
</script>

<HeatmapTable
  data={metrics_data}
  columns={show_cols}
  higher_is_better={geo_opt.higher_is_better}
  lower_is_better={geo_opt.lower_is_better}
  format={{ RMSD: `.3f`, Match: `.1%`, Decrease: `.1%`, Increase: `.1%` }}
/>

<style>
  :global(.heatmap-table td:not(:first-child)) {
    text-align: right;
  }
</style>
