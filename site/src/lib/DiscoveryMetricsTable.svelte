<script lang="ts">
  import { HeatmapTable, MODEL_METADATA, TRAINING_SETS, model_is_compliant } from '$lib'
  import { discovery, phonons } from '$root/scripts/metrics-which-is-better.yml'
  import { pretty_num } from 'elementari'

  export let discovery_set: `full_test_set` | `most_stable_10k` | `unique_prototypes` =
    `unique_prototypes`
  export let show_non_compliant: boolean = false
  export let show_energy_only: boolean = false
  export let show_metadata: boolean = true
  export let metadata_cols = [
    { label: `Training Set` },
    { label: `Params` },
    { label: `Targets` },
    { label: `Date Added` },
  ]
  export let hide_cols: string[] = []
  export let columns = [
    { label: `Model` },
    { label: `F1`, tooltip: `harmonic mean of precision and recall` },
    { label: `DAF`, tooltip: `discovery acceleration factor` },
    { label: `Prec`, tooltip: `precision of classifying thermodynamic stability` },
    { label: `Acc`, tooltip: `accuracy of classifying thermodynamic stability` },
    {
      label: `TPR`,
      tooltip: `true positive rate of classifying thermodynamic stability`,
    },
    {
      label: `TNR`,
      tooltip: `true negative rate of classifying thermodynamic stability`,
    },
    {
      label: `MAE`,
      tooltip: `mean absolute error of predicting the convex hull distance`,
    },
    { label: `RMSE` },
    { label: `R<sup>2</sup>`, tooltip: `coefficient of determination` },
    {
      label: `κ<sub>SRME</sub>`,
      tooltip: `symmetric relative mean error in predicted phonon mode contributions to thermal conductivity κ`,
    },
    ...(show_metadata ? metadata_cols : []),
  ].filter((col) => !hide_cols.includes(col.label))
  export let sep_lines: number[] = [6, 9]

  function format_train_set(model_training_sets: string[]) {
    let [total_structs, total_materials] = [0, 0]
    const data_urls: Record<string, string> = {}
    const tooltip: string[] = []

    for (const train_set of model_training_sets) {
      if (!(train_set in TRAINING_SETS)) {
        console.warn(`Training set ${train_set} not found in TRAINING_SETS`)
        continue
      }
      const training_set_info = TRAINING_SETS[train_set]
      const n_structs = training_set_info.n_structures
      const n_materials = training_set_info.n_materials ?? n_structs

      total_structs += n_structs
      total_materials += n_materials

      const title = training_set_info.title || train_set
      data_urls[train_set || title] = training_set_info.url || ``

      if (n_materials !== n_structs) {
        tooltip.push(
          `${title}: ${pretty_num(n_materials, `,`)} materials (${pretty_num(n_structs, `,`)} structures)`,
        )
      } else {
        tooltip.push(`${title}: ${pretty_num(n_materials, `,`)} materials`)
      }
    }

    const data_links = Object.entries(data_urls).map(([key, href]) => {
      if (href) {
        return `<a href="${href}" target="_blank" rel="noopener noreferrer">${key}</a>`
      }
      return key
    })

    const data_str = data_links.join(`+`)
    const new_line = `&#013;` // line break that works in title attribute
    const dataset_tooltip =
      tooltip.length > 1 ? `${new_line}• ${tooltip.join(new_line + `• `)}` : ``

    let title = `${pretty_num(total_materials, `,`)} materials in training set${new_line}${dataset_tooltip}`
    let train_size_str = `<span title="${title}" data-sort-value="${total_materials}">${pretty_num(total_materials)} (${data_str})</span>`

    if (total_materials !== total_structs) {
      title =
        `${pretty_num(total_materials, `,`)} materials in training set ` +
        `(${pretty_num(total_structs, `,`)} structures counting all DFT relaxation ` +
        `frames per material)${dataset_tooltip}`

      train_size_str =
        `<span title="${title}" data-sort-value="${total_materials}">` +
        `${pretty_num(total_materials)} <small>(${pretty_num(total_structs)})</small> ` +
        `(${data_str})</span>`
    }

    return train_size_str
  }

  const Targets = {
    E: `Energy`,
    EF_C: `Energy with conservative forces`,
    EF_D: `Energy with direct forces`,
    EFS_C: `Energy with conservative forces and stress`,
    EFS_D: `Energy with direct forces and stress`,
    EFS_CM: `Energy with conservative forces, stress, and Magmoms`,
    EFS_DM: `Energy with direct forces, stress, and magmoms`,
  }

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
      (show_energy_only || model.targets != `E`) &&
      (show_non_compliant || model_is_compliant(model)),
  )
    .map((model) => {
      const disc_metrics = model.metrics?.discovery?.[discovery_set]

      const targets = model.targets.replace(/_(.)/g, `<sub>$1</sub>`)
      const targets_str = `<span title="${Targets[model.targets]}">${targets}</span>`

      // rename metric keys to pretty labels
      return {
        Model: `<a title="Version: ${model.model_version}" href="/models/${model.model_key}">${model.model_name}</a>`,
        F1: disc_metrics?.F1,
        DAF: disc_metrics?.DAF,
        Prec: disc_metrics?.Precision,
        Acc: disc_metrics?.Accuracy,
        TPR: disc_metrics?.TPR,
        TNR: disc_metrics?.TNR,
        MAE: disc_metrics?.MAE,
        RMSE: disc_metrics?.RMSE,
        'R<sup>2</sup>': disc_metrics?.R2,
        'κ<sub>SRME</sub>': model.metrics?.phonons?.κ_SRME,
        'Training Set': format_train_set(model.training_set),
        Params: `<span title="${pretty_num(model.model_params, `,`)} trainable model parameters">${pretty_num(model.model_params)}</span>`,
        Targets: targets_str,
        'Date Added': `<span title="${long_date(model.date_added)}">${model.date_added}</span>`,
      }
    })
    .sort((row1, row2) => (row2.F1 ?? 0) - (row1.F1 ?? 0)) // Sort by F1 score descending
</script>

<HeatmapTable
  data={metrics_data}
  {columns}
  higher_is_better={[...discovery.higher_is_better, ...phonons.higher_is_better]}
  lower_is_better={[...discovery.lower_is_better, ...phonons.lower_is_better]}
  {sep_lines}
/>
