<script lang="ts">
  import { HeatmapTable, MODEL_METADATA, model_is_compliant } from '$lib'
  import { geo_opt } from '$root/scripts/metrics-which-is-better.yml'

  export let show_non_compliant: boolean = false
  export let show_metadata: boolean = true
  export let metadata_cols: { label: string; tooltip?: string }[] = []
  export let columns: { label: string; tooltip?: string; style?: string }[] = [
    { label: `Model` },
    {
      label: `RMSD`,
      tooltip: `Root mean squared displacement (in Å) of ML vs DFT relaxed atomic positions as calculated by pymatgen StructureMatcher`,
    },
    {
      label: `σ<sub>match</sub>`,
      tooltip: `Fraction of structures where ML and DFT ground state have matching spacegroup`,
    },
    {
      label: `σ<sub>dec</sub>`,
      tooltip: `Fraction of structures where the number of symmetry operations decreased after ML relaxation`,
    },
    {
      label: `σ<sub>inc</sub>`,
      tooltip: `Fraction of structures where the number of symmetry operations increased after ML relaxation`,
    },
    // { label: `n_structs` },
    ...(show_metadata ? metadata_cols : []),
  ]

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
        'σ<sub>match</sub>': geo_opt.symmetry_match,
        'σ<sub>dec</sub>': geo_opt.symmetry_decrease,
        'σ<sub>inc</sub>': geo_opt.symmetry_increase,
        // n_structs: geo_opt.n_structs,
        'Date Added': `<span title="${long_date(model.date_added)}">${model.date_added}</span>`,
      }
    })
    // Sort by Match descending
    .sort(
      (row1, row2) =>
        (row2?.[`σ<sub>match</sub>`] ?? 0) - (row1?.[`σ<sub>match</sub>`] ?? 0),
    )
</script>

<HeatmapTable
  data={metrics_data}
  {columns}
  higher_is_better={geo_opt.higher_is_better}
  lower_is_better={geo_opt.lower_is_better}
  format={{
    RMSD: `.3f`,
    'σ<sub>match</sub>': `.1%`,
    'σ<sub>dec</sub>': `.1%`,
    'σ<sub>inc</sub>': `.1%`,
  }}
  {...$$restProps}
/>

<style>
  :global(.heatmap-table td:not(:first-child)) {
    text-align: right;
  }
</style>
