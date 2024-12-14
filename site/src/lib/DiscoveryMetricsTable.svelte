<script lang="ts">
  import {
    HeatmapTable,
    MODEL_METADATA,
    TRAINING_SETS,
    get_metric_rank_order,
    model_is_compliant,
  } from '$lib'
  import { pretty_num } from 'elementari'
  import { Tooltip } from 'svelte-zoo'
  import type { HeatmapColumn, ModelData } from './types.ts'

  export let discovery_set: `full_test_set` | `most_stable_10k` | `unique_prototypes` =
    `unique_prototypes`
  export let show_non_compliant: boolean = false
  export let show_energy_only: boolean = false
  export let show_metadata: boolean = true
  export let hide_cols: string[] = []
  export let metadata_cols = [
    { label: `Training Set`, tooltip: `Size of and link to model training set` },
    { label: `Params`, tooltip: `Number of trainable model parameters` },
    { label: `Targets`, tooltip: `Target property used to train the model` },
    { label: `Date Added`, tooltip: `Submission date to the leaderboard` },
  ]

  export let f1_weight = 0.5
  $: kappa_weight = 1 - f1_weight

  let columns: HeatmapColumn[]
  $: columns = [
    { label: `Model`, sticky: true },
    {
      label: `E<sub>MBD</sub>`,
      tooltip: `Combined effectiveness score: ${(f1_weight * 100).toFixed()}% F1 + ${(kappa_weight * 100).toFixed()}% κ<sub>SRME</sub>`,
      style: `border-left: 1px solid black;`,
      format: `.3f`,
      better: `higher`,
    },
    { label: `F1`, tooltip: `Harmonic mean of precision and recall` },
    { label: `DAF`, tooltip: `Discovery acceleration factor` },
    { label: `Prec`, tooltip: `Precision of classifying thermodynamic stability` },
    { label: `Acc`, tooltip: `Accuracy of classifying thermodynamic stability` },
    {
      label: `TPR`,
      tooltip: `True positive rate of classifying thermodynamic stability`,
    },
    {
      label: `TNR`,
      tooltip: `True negative rate of classifying thermodynamic stability`,
    },
    {
      label: `MAE`,
      tooltip: `Mean absolute error of predicting the convex hull distance`,
      style: `border-left: 1px solid black;`,
    },
    {
      label: `RMSE`,
      tooltip: `Root mean squared error of predicting the convex hull distance`,
    },
    { label: `R<sup>2</sup>`, tooltip: `Coefficient of determination` },
    {
      label: `κ<sub>SRME</sub>`,
      tooltip: `Symmetric relative mean error in predicted phonon mode contributions to thermal conductivity κ`,
      style: `border-left: 1px solid black;`,
    },
    ...(show_metadata ? metadata_cols : []),
  ].map((col) => ({
    ...col,
    better: col.better ?? get_metric_rank_order(col.label),
    hidden: hide_cols.includes(col.label) || col.hidden,
  }))

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

  const targets_tooltips: Record<ModelData[`targets`], string> = {
    E: `Energy`,
    EF_G: `Energy with gradient-based forces`,
    EF_D: `Energy with direct forces`,
    EFS_G: `Energy with gradient-based forces and stress`,
    EFS_D: `Energy with direct forces and stress`,
    EFS_GM: `Energy with gradient-based forces, stress, and magmoms`,
    EFS_DM: `Energy with direct forces, stress, and magmoms`,
  }

  const long_date = (date: string): string =>
    new Date(date).toLocaleDateString(undefined, {
      weekday: `long`,
      year: `numeric`,
      month: `long`,
      day: `numeric`,
    })

  // Constants for test set sizes and derived weights
  const N_WBM = 257000
  const N_PHONONS = 103
  const WBM_WEIGHT = Math.log(N_WBM)
  const PHONONS_WEIGHT = Math.log(N_PHONONS)

  // Beta for F-score calculation (β = 1/2 weights precision twice as much as recall)
  const BETA = 1 / 2
  const BETA_SQ = BETA * BETA

  // Calculate F1 score from precision and recall
  function calculate_f_beta(precision: number, recall: number): number {
    return ((1 + BETA_SQ) * (precision * recall)) / (BETA_SQ * precision + recall)
  }

  // Convert error metric (SRME) to effectiveness (1 - normalized_error)
  function error_to_effectiveness(error: number): number {
    return 1 - Math.min(1, error) // Clamp to [0,1] range
  }

  // Calculate combined effectiveness score
  $: calculate_effectiveness = (metrics: any, kappa: number | undefined) => {
    if (!metrics?.Precision || !metrics?.Recall || !kappa) return undefined

    const f_beta = calculate_f_beta(metrics.Precision, metrics.Recall)
    const kappa_effectiveness = error_to_effectiveness(kappa)

    // Use f1_weight from slider for weighted average
    return f1_weight * f_beta + kappa_weight * kappa_effectiveness
  }

  // Transform MODEL_METADATA into table data format
  $: metrics_data = MODEL_METADATA.filter(
    (model) =>
      (show_energy_only || model.targets != `E`) &&
      (show_non_compliant || model_is_compliant(model)),
  )
    .map((model) => {
      const metrics = model.metrics?.discovery?.[discovery_set]
      const kappa = model.metrics?.phonons?.κ_SRME

      // Calculate F1 score
      const f_beta =
        metrics?.Precision && metrics?.Recall
          ? calculate_f_beta(metrics.Precision, metrics.Recall)
          : undefined

      // Calculate overall effectiveness score
      const effectiveness = calculate_effectiveness(metrics, kappa)

      return {
        Model: `<a title="Version: ${model.model_version}" href="/models/${model.model_key}">${model.model_name}</a>`,
        'E<sub>MBD</sub>': effectiveness,
        F1: metrics?.F1,
        DAF: metrics?.DAF,
        Prec: metrics?.Precision,
        Acc: metrics?.Accuracy,
        TPR: metrics?.TPR,
        TNR: metrics?.TNR,
        MAE: metrics?.MAE,
        RMSE: metrics?.RMSE,
        'R<sup>2</sup>': metrics?.R2,
        'κ<sub>SRME</sub>': kappa,
        'Training Set': format_train_set(model.training_set),
        Params: `<span title="${pretty_num(model.model_params, `,`)} trainable model parameters">${pretty_num(model.model_params)}</span>`,
        Targets: `<span title="${targets_tooltips[model.targets]}">${model.targets.replace(/_(.)/g, `<sub>$1</sub>`)}</span>`,
        'Date Added': `<span title="${long_date(model.date_added)}">${model.date_added}</span>`,
      }
    })
    .sort((row1, row2) => (row2[`E<sub>MBD</sub>`] ?? 0) - (row1[`E<sub>MBD</sub>`] ?? 0))

  // Update tooltip to show current weights
  $: columns = columns?.map((col) =>
    col.label === `E<sub>MBD</sub>`
      ? {
          ...col,
          tooltip: `Combined effectiveness score: ${(f1_weight * 100).toFixed()}% F1 + ${(kappa_weight * 100).toFixed()}% κ<sub>SRME</sub>`,
        }
      : col,
  )
</script>

<div class="weight-slider">
  <Tooltip>
    <span slot="tip">Adjust relative weights of F1 and κ<sub>SRME</sub></span>
    <label>
      F1 {(f1_weight * 100).toFixed()}% - κ<sub>SRME</sub>
      {(kappa_weight * 100).toFixed()}%
      <input type="range" min="0" max="1" step="0.1" bind:value={f1_weight} />
    </label>
  </Tooltip>
</div>

<HeatmapTable data={metrics_data} {columns} {...$$restProps} />

<style>
  .weight-slider {
    display: flex;
    justify-content: center;
    margin-bottom: 1em;
  }
  .weight-slider label {
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 0.5em;
    font-size: 0.9em;
  }
  .weight-slider input {
    width: 200px;
  }
</style>
