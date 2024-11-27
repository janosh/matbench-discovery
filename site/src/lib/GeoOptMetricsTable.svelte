<script lang="ts">
  import { HeatmapTable, MODEL_METADATA, model_is_compliant } from '$lib'
  import { geo_opt } from '$root/scripts/metrics-which-is-better.yml'
  import { pretty_num } from 'elementari'

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
    {
      label: `N<sub>structs</sub>`,
      tooltip: `Number of structures relaxed by each model and used for these metrics`,
    },
    ...(show_metadata ? metadata_cols : []),
  ]

  // Transform MODEL_METADATA into table data format
  $: metrics_data = MODEL_METADATA.filter(
    (model) =>
      (show_non_compliant || model_is_compliant(model)) &&
      model.metrics?.geo_opt?.rmsd != undefined,
  )
    .sort(
      (row1, row2) =>
        (row2?.metrics?.geo_opt?.symmetry_match ?? 0) -
        (row1?.metrics?.geo_opt?.symmetry_match ?? 0),
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
        'N<sub>structs</sub>': pretty_num(geo_opt.n_structures),
      }
    })
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
