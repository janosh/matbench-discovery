<script lang="ts">
  import {
    HeatmapTable,
    MODEL_METADATA,
    model_is_compliant,
    get_metric_rank_order,
  } from '$lib'
  import { pretty_num } from 'elementari'
  import type { HeatmapColumn } from './types.ts'

  export let show_non_compliant: boolean = false
  export let show_metadata: boolean = true
  export let metadata_cols: { label: string; tooltip?: string }[] = []

  // Get all unique symprec values from MODEL_METADATA
  $: symprec_values = [
    ...new Set(
      MODEL_METADATA.flatMap((model) =>
        Object.keys(model.metrics?.geo_opt ?? {})
          .filter((key) => key.startsWith(`symprec=`))
          .map((key) => key.replace(`symprec=`, ``)),
      ),
    ),
  ].sort((val1, val2) => parseFloat(val2) - parseFloat(val1)) // Sort in descending order

  // Helper to format symprec in scientific notation
  const format_symprec = (symprec: string) => {
    const exp = symprec.split(`e-`)[1]
    return `symprec=10<sup>-${exp}</sup>Å`
  }

  const sep_line_style = `border-left: 1px solid black`

  // Create columns for each symprec value
  let columns: HeatmapColumn[]
  $: columns = [
    { label: `Model`, sticky: true },
    {
      label: `RMSD`,
      tooltip: `Root mean squared displacement (in Å) of ML vs DFT relaxed atomic positions as calculated by pymatgen StructureMatcher`,
      style: sep_line_style,
      format: `.3f`,
    },
    // Symmetry match columns
    ...symprec_values.flatMap((symprec) => [
      {
        group: format_symprec(symprec),
        label: `σ<sub>match</sub>`,
        tooltip: `Fraction of structures where ML and DFT ground state have matching spacegroup (symprec=${symprec}Å)`,
        style: sep_line_style,
        format: `.1%`,
      },
      {
        group: format_symprec(symprec),
        label: `σ<sub>dec</sub>`,
        tooltip: `Fraction of structures where the number of symmetry operations decreased after ML relaxation (symprec=${symprec}Å)`,
        format: `.1%`,
      },
      {
        group: format_symprec(symprec),
        label: `σ<sub>inc</sub>`,
        tooltip: `Fraction of structures where the number of symmetry operations increased after ML relaxation (symprec=${symprec}Å). Not colored because it's high or low is good or bad. Could be models find higher symmetry lower-energy structures than DFT optimizer.`,
        color_scale: null,
        format: `.1%`,
      },
      {
        group: format_symprec(symprec),
        label: `N<sub>ops,MAE</sub>`,
        tooltip: `Mean absolute error of number of symmetry operations in DFT and ML-relaxed structures (symprec=${symprec}Å)`,
        format: `.3`,
      },
    ]),
    {
      label: `N<sub>structs</sub>`,
      tooltip: `Number of structures relaxed by each model and used to compute these metrics`,
      style: sep_line_style,
    },
    ...(show_metadata ? metadata_cols : []),
  ].map((col) => ({ ...col, better: col.better ?? get_metric_rank_order(col.label) }))

  // Transform MODEL_METADATA into table data format
  $: metrics_data = MODEL_METADATA.filter(
    (model) =>
      (show_non_compliant || model_is_compliant(model)) &&
      // Check if model has data for all symprec values
      symprec_values.every(
        (symprec) => model.metrics?.geo_opt?.[`symprec=${symprec}`]?.rmsd != undefined,
      ) &&
      model.model_name !== `BOWSR`, // hide BOWSR as it's a huge outlier that makes the table hard to read
  )
    .sort(
      (row1, row2) =>
        (row2?.metrics?.geo_opt?.[`symprec=${symprec_values[0]}`]?.symmetry_match ?? 0) -
        (row1?.metrics?.geo_opt?.[`symprec=${symprec_values[0]}`]?.symmetry_match ?? 0),
    )
    .map((model) => {
      const geo_opt = model.metrics?.geo_opt
      if (!geo_opt) return null

      return {
        Model: `<a title="Version: ${model.model_version}" href="/models/${model.model_key}">${model.model_name}</a>`,
        RMSD: geo_opt[`symprec=${symprec_values[0]}`].rmsd,
        ...symprec_values.reduce(
          (acc, symprec) => ({
            ...acc,
            [`σ<sub>match</sub> (${format_symprec(symprec)})`]:
              geo_opt[`symprec=${symprec}`].symmetry_match,
            [`σ<sub>dec</sub> (${format_symprec(symprec)})`]:
              geo_opt[`symprec=${symprec}`].symmetry_decrease,
            [`σ<sub>inc</sub> (${format_symprec(symprec)})`]:
              geo_opt[`symprec=${symprec}`].symmetry_increase,
            [`N<sub>ops,MAE</sub> (${format_symprec(symprec)})`]:
              geo_opt[`symprec=${symprec}`].n_sym_ops_mae,
          }),
          {},
        ),
        'N<sub>structs</sub>': `<span title="${model.model_name} relaxed ${pretty_num(
          geo_opt[`symprec=${symprec_values[0]}`].n_structures,
          `,`,
        )} structures">${pretty_num(geo_opt[`symprec=${symprec_values[0]}`].n_structures)}</span>`,
      }
    })
</script>

<HeatmapTable data={metrics_data} {columns} {...$$restProps} />

<style>
  :global(.heatmap-table td:not(:first-child)) {
    text-align: right;
  }
</style>
