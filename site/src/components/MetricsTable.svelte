<script lang="ts">
  import { MODEL_METADATA, TRAINING_SETS } from '$lib'
  import { si_fmt } from '$lib/utils'
  import { pretty_num } from 'elementari'
  import HeatmapTable from './HeatmapTable.svelte'

  export let discovery_set: `full` | `most_stable_10k` | `unique_prototypes` =
    `unique_prototypes`

  const higherIsBetter = new Set([
    `F1`,
    `DAF`,
    `Prec`,
    `Acc`,
    `TPR`,
    `TNR`,
    `R<sup>2</sup>`,
  ])

  const lowerIsBetter = new Set([`MAE`, `RMSE`, `κ<sub>SRME</sub>`])

  const metaColumns = [`Training Set`, `Params`, `Model Type`, `Targets`, `Date Added`]

  const show_cols = [
    `Model`,
    `F1`,
    `DAF`,
    `Prec`,
    `Acc`,
    `TPR`,
    `TNR`,
    `MAE`,
    `RMSE`,
    `R<sup>2</sup>`,
    `κ<sub>SRME</sub>`,
    ...(discovery_set === `unique_prototypes` ? metaColumns : []),
  ]

  function format_train_set(model_training_sets: string[]) {
    let [total_structs, total_materials] = [0, 0]
    const data_urls: Record<string, string> = {}
    const tooltipLines: string[] = []

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
        tooltipLines.push(
          `${title}: ${si_fmt(n_materials)} materials (${si_fmt(n_structs)} structures)`,
        )
      } else {
        tooltipLines.push(`${title}: ${si_fmt(n_materials)} materials`)
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
      tooltipLines.length > 1 ? `${new_line}• ${tooltipLines.join(new_line + `• `)}` : ``

    let title = `${pretty_num(total_materials)} materials in training set${new_line}${dataset_tooltip}`
    let trainSizeStr = `<span title="${title}" data-sort-value="${total_materials}">${si_fmt(total_materials)} (${data_str})</span>`

    if (total_materials !== total_structs) {
      title =
        `${pretty_num(total_materials)} materials in training set ` +
        `(${pretty_num(total_structs)} structures counting all DFT relaxation ` +
        `frames per material)${dataset_tooltip}`

      trainSizeStr =
        `<span title="${title}" data-sort-value="${total_materials}">` +
        `${si_fmt(total_materials)} <small>(${si_fmt(total_structs)})</small> ` +
        `(${data_str})</span>`
    }

    return trainSizeStr
  }

  const long_date = (date: string): string =>
    new Date(date).toLocaleDateString(undefined, {
      weekday: `long`,
      year: `numeric`,
      month: `long`,
      day: `numeric`,
    })

  // Transform MODEL_METADATA into table data format
  $: metrics_data = MODEL_METADATA.map((model) => {
    const disc_metrics = model.metrics?.discovery?.[discovery_set]

    // rename metric keys to pretty labels
    return {
      Model: `<span title="Version: ${model.model_version}">${model.model_name}</span>`,
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
      Params: `<span title="${pretty_num(model.model_params, `,`)}">${pretty_num(model.model_params)}</span>`,
      Targets: model.targets.replace(/_(.)/g, `<sub>$1</sub>`),
      'Date Added': `<span title="${long_date(model.date_added)}">${model.date_added}</span>`,
    }
  }).sort((a, b) => (b.F1 ?? 0) - (a.F1 ?? 0)) // Sort by F1 score descending
</script>

<HeatmapTable data={metrics_data} columns={show_cols} {higherIsBetter} {lowerIsBetter} />
