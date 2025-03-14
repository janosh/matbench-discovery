<script lang="ts">
  import { HeatmapTable, MODEL_METADATA, model_is_compliant } from '$lib'
  import { pretty_num } from 'elementari'
  import type { HeatmapColumn } from './types.ts'

  interface Props {
    show_non_compliant?: boolean
  }

  let { show_non_compliant = false, ...rest }: Props = $props()

  // Get all unique symprec values from MODEL_METADATA
  let symprec_values = $derived(
    [
      ...new Set(
        MODEL_METADATA.flatMap((model) =>
          Object.keys(model.metrics?.geo_opt ?? {})
            .filter((key) => key.startsWith(`symprec=`))
            .map((key) => key.replace(`symprec=`, ``)),
        ),
      ),
    ].sort((val1, val2) => parseFloat(val2) - parseFloat(val1)),
  ) // Sort in descending order

  // Helper to format symprec in scientific notation
  const format_symprec = (symprec: string) => `10<sup>-${symprec.split(`e-`)[1]}</sup>Å`

  const sep_line_style = `border-left: 1px solid black`

  // Define base columns
  const base_cols: HeatmapColumn[] = [
    { label: `Model`, sticky: true },
    {
      label: `RMSD`,
      tooltip: `Root mean squared displacement (in Å) of ML vs DFT relaxed atomic positions as calculated by pymatgen StructureMatcher`,
      style: sep_line_style,
      format: `.3f`,
    },
  ]

  // Define symmetry metrics
  const sym_metrics = [
    {
      key: `match`,
      label: `σ<sub>match</sub>`,
      tooltip: (symprec: string) =>
        `Fraction of structures where ML and DFT ground state have matching spacegroup at ${format_symprec(symprec)}`,
    },
    {
      key: `dec`,
      label: `σ<sub>dec</sub>`,
      tooltip: (symprec: string) =>
        `Fraction of structures where the number of symmetry operations decreased after ML relaxation at ${format_symprec(symprec)}`,
    },
    {
      key: `inc`,
      label: `σ<sub>inc</sub>`,
      tooltip: (symprec: string) =>
        `Fraction of structures where the number of symmetry operations increased after ML relaxation at ${format_symprec(symprec)}. Not colored as high/low values are neither good nor bad.`,
      color_scale: undefined,
    },
    {
      key: `ops_mae`,
      label: `N<sub>ops,MAE</sub>`,
      tooltip: (symprec: string) =>
        `Mean absolute error of number of symmetry operations in DFT and ML-relaxed structures at ${format_symprec(symprec)}`,
    },
  ]

  // Create columns for each symprec value
  let columns = $derived([
    ...base_cols,
    ...symprec_values.flatMap((symprec) =>
      sym_metrics.map(
        ({ key, label, tooltip, ...rest }, idx): HeatmapColumn => ({
          group: format_symprec(symprec),
          label,
          tooltip: tooltip(symprec),
          style: idx === 0 ? sep_line_style : undefined,
          format: key === `ops_mae` ? `.3` : `.1%`,
          ...rest,
        }),
      ),
    ),
    {
      label: `N<sub>structs</sub>`,
      tooltip: `Number of structures relaxed by each model and used to compute these metrics`,
      style: sep_line_style,
    },
    {
      label: `Links`,
      tooltip: `Links to model resources`,
      style: sep_line_style,
      format: undefined,
      better: undefined,
    },
  ])

  // Transform MODEL_METADATA into table data format
  let metrics_data = $derived(
    MODEL_METADATA.filter(
      (model) =>
        (show_non_compliant || model_is_compliant(model)) &&
        // Check if model has data for all symprec values
        symprec_values.every(
          (symprec) => model.metrics?.geo_opt?.[`symprec=${symprec}`]?.rmsd != undefined,
        ) &&
        model.model_name !== `BOWSR`, // hide BOWSR as it's a huge outlier that makes the table hard to read
    )
      .map((model) => {
        const geo_opt = model.metrics?.geo_opt
        if (!geo_opt) return undefined

        const symprec_key = `symprec=${symprec_values[0]}`
        const result = {
          Model: `<a title="Version: ${model.model_version}" href="/models/${model.model_key}">${model.model_name}</a>`,
          RMSD: geo_opt[symprec_key]?.rmsd,
          ...symprec_values.reduce((acc, symprec) => {
            const metrics = geo_opt[`symprec=${symprec}`]
            if (!metrics) return acc

            return {
              ...acc,
              [`σ<sub>match</sub> (${format_symprec(symprec)})`]: metrics.symmetry_match,
              [`σ<sub>dec</sub> (${format_symprec(symprec)})`]: metrics.symmetry_decrease,
              [`σ<sub>inc</sub> (${format_symprec(symprec)})`]: metrics.symmetry_increase,
              [`N<sub>ops,MAE</sub> (${format_symprec(symprec)})`]: metrics.n_sym_ops_mae,
            }
          }, {}),
          'N<sub>structs</sub>': `<span title="${model.model_name} relaxed ${pretty_num(
            geo_opt[symprec_key]?.n_structures ?? 0,
            `,`,
          )} structures">${pretty_num(geo_opt[symprec_key]?.n_structures ?? 0)}</span>`,
          Links: {
            files: [
              {
                url: geo_opt.pred_file_url,
                title: `Download model-relaxed WBM structures`,
                icon: `📦`,
              },
              ...symprec_values.map((symprec) => ({
                url: geo_opt[`symprec=${symprec}`]?.analysis_file_url,
                title: `Download ${model.model_name}-relaxed WBM structure analysis for symprec ${format_symprec(symprec)}`,
                icon: `📊<sup>-${symprec.split(`e-`)[1]}</sup>`,
              })),
            ],
          },
        }
        return result
      })
      .filter((row) => row !== undefined),
  )
</script>

<HeatmapTable data={metrics_data} {columns} {...rest}>
  {#snippet cell({ col, val })}
    {#if col.label === `Links` && val}
      {@const links = val}
      {#each links.files as { url: href, title, icon } (href)}
        {#if href}
          <a {href} {title} target="_blank" rel="noopener noreferrer">
            {@html icon}
          </a>
        {/if}
      {/each}
    {:else if typeof val === `number` && col.format}
      {pretty_num(val, col.format)}
    {:else if [undefined, null].includes(val)}
      n/a
    {:else}
      {@html val}
    {/if}
  {/snippet}
</HeatmapTable>
