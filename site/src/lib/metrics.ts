import {
  DATASETS,
  MODELS,
  format_date,
  get_pred_file_urls,
  model_is_compliant,
} from '$lib'
import type { TargetType } from '$lib/model-schema'
import modeling_tasks from '$pkg/modeling-tasks.yml'
import { max, min } from 'd3-array'
import { scaleLog, scaleSequential } from 'd3-scale'
import * as d3sc from 'd3-scale-chromatic'
import { choose_bw_for_contrast, pretty_num } from 'elementari/labels'
import { DEFAULT_CPS_CONFIG, calculate_cps } from './combined_perf_score'
import type { DiscoveryMetricsSet } from './model-schema.d.ts'
import type {
  CombinedMetricConfig,
  DiscoverySet,
  HeatmapColumn,
  LinkData,
  ModelData,
  ModelStats,
  RowData,
} from './types'

// model target type descriptions
export const targets_tooltips: Record<TargetType, string> = {
  E: `Energy`,
  EF_G: `Energy with gradient-based forces`,
  EF_D: `Energy with direct forces`,
  EFS_G: `Energy with gradient-based forces and stress`,
  EFS_D: `Energy with direct forces and stress`,
  EFS_GM: `Energy with gradient-based forces, stress, and magmoms`,
  EFS_DM: `Energy with direct forces, stress, and magmoms`,
}

// Gets a metric value from a model, handling metrics from different tasks
// and different locations within the model data structure
export function get_metric_value(
  model: ModelData, // model data to extract the metric from
  metric_key: string, // key of the metric to extract
  set: string = `full_test_set`, // optional dataset to use
): number | undefined {
  // Special case for CPS
  if (metric_key === `CPS`) {
    return model.CPS
  }

  // Check if it's a discovery metric
  if (model.metrics?.discovery?.[set as keyof typeof model.metrics.discovery]) {
    const metrics_set = model.metrics.discovery[
      set as keyof typeof model.metrics.discovery
    ] as DiscoveryMetricsSet
    if (
      metrics_set &&
      metrics_set[metric_key as keyof DiscoveryMetricsSet] !== undefined
    ) {
      const value = metrics_set[metric_key as keyof DiscoveryMetricsSet]
      return typeof value === `number` ? value : Number(value)
    }
  }

  if (
    [`κ_SRME`, `kappa_SRME`].includes(metric_key) &&
    model.metrics?.phonons &&
    typeof model.metrics.phonons !== `string`
  ) {
    const value = model.metrics.phonons.kappa_103?.κ_SRME
    return value !== undefined ? Number(value) : undefined
  }

  // Check for RMSD in geo_opt
  if (
    metric_key === `RMSD` &&
    model.metrics?.geo_opt &&
    typeof model.metrics.geo_opt !== `string`
  ) {
    return model.metrics.geo_opt[`symprec=1e-5`]?.rmsd
  }

  // Try to access any other potential location using dotted path notation
  const nested_val = metric_key.split(`.`).reduce<unknown>((acc, key) => {
    if (acc && typeof acc === `object`) {
      return (acc as Record<string, unknown>)[key]
    }
    return undefined
  }, model.metrics)

  if (nested_val !== undefined) {
    return typeof nested_val === `number` ? nested_val : Number(nested_val)
  }

  return undefined
}

export const all_higher_better_metrics = Object.values(modeling_tasks).flatMap(
  (model_task) => model_task.metrics.higher_is_better,
)
export const all_lower_better_metrics = Object.values(modeling_tasks).flatMap(
  (model_task) => model_task.metrics.lower_is_better,
)
export function metric_better_as(metric: string) {
  if (all_higher_better_metrics.includes(metric)) return `higher`
  if (all_lower_better_metrics.includes(metric)) return `lower`
  return null
}

// Format training set information for display in the metrics table
export function format_train_set(model_train_sets: string[], model: ModelData): string {
  const { n_training_structures = 0, n_training_materials = 0 } = model

  const data_urls: Record<string, string> = {}
  const tooltip: string[] = []

  for (const data_name of model_train_sets) {
    if (!(data_name in DATASETS)) {
      console.warn(`Training set ${data_name} not found in DATASETS`)
      continue
    }
    const { title, slug, n_structures, n_materials = n_structures } = DATASETS[data_name]
    data_urls[data_name] = `/data/${slug}`

    if (n_materials !== n_structures) {
      tooltip.push(
        `${title}: ${pretty_num(n_materials, `,`)} materials (${pretty_num(n_structures, `,`)} structures)`,
      )
    } else {
      tooltip.push(`${title}: ${pretty_num(n_materials, `,`)} materials`)
    }
  }

  const dataset_links = Object.entries(data_urls)
    .map(([key, href]) => `<a href="${href}">${key}</a>`)
    .join(`+`)
  const new_line = `&#013;` // line break that works in title attribute
  const dataset_tooltip =
    tooltip.length > 1 ? `${new_line}• ${tooltip.join(new_line + `• `)}` : ``

  let title = `${pretty_num(n_training_materials, `,`)} materials in training set${new_line}${dataset_tooltip}`
  let train_size_str = `<span title="${title}" data-sort-value="${n_training_materials}">${pretty_num(n_training_materials)} <small>${dataset_links}</small></span>`

  if (n_training_materials !== n_training_structures) {
    title =
      `${pretty_num(n_training_materials, `,`)} materials in training set ` +
      `(${pretty_num(n_training_structures, `,`)} structures counting all DFT relaxation ` +
      `frames per material)${dataset_tooltip}`

    train_size_str =
      `<span title="${title}" data-sort-value="${n_training_materials || n_training_structures}">` +
      `${pretty_num(n_training_materials)} <small>(${pretty_num(n_training_structures)})</small> ` +
      `<small>${dataset_links}</small></span>`
  }

  return train_size_str
}

// Safely access nested geometry optimization properties
export function get_geo_opt_property<T>(
  geo_opt: unknown, // geometry optimization data
  symprec: string, // symmetry precision value
  property: string, // property name to retrieve
): T | undefined {
  if (!geo_opt || typeof geo_opt !== `object`) return undefined

  const symprec_key = `symprec=${symprec}`
  const metrics = geo_opt[symprec_key as keyof typeof geo_opt]

  if (!metrics || typeof metrics !== `object`) return undefined

  return metrics[property as keyof typeof metrics] as T
}

// combined filter function that respects both prediction type and compliance filters
export function make_combined_filter(
  model_filter: (model: ModelData) => boolean, // custom model filter function
  show_energy: boolean, // whether to show energy-only models
  show_noncomp: boolean, // whether to show non-compliant models
): (model: ModelData) => boolean {
  return (model: ModelData) => {
    if (!model_filter(model)) return false // Apply user-provided model_filter first

    // Filter energy-only models if not shown
    const is_energy_only = model.targets === `E`
    if (is_energy_only && !show_energy) return false

    // Filter noncompliant models if not shown
    const is_compliant = model_is_compliant(model)
    if (!is_compliant && !show_noncomp) return false

    return true
  }
}

// Calculate table data for the metrics table with combined scores
export function calculate_metrics_data(
  discovery_set: DiscoverySet, // discovery set to use for metrics
  model_filter: (model: ModelData) => boolean, // filter function for models
  show_energy_only: boolean, // show energy-only models
  show_noncompliant: boolean, // show non-compliant models
  config: CombinedMetricConfig, // combined metric configuration
  compliant_clr: string = `#4caf50`, // color for compliant models
  noncompliant_clr: string = `#4682b4`, // color for non-compliant models
): RowData[] {
  const current_filter = make_combined_filter(
    model_filter,
    show_energy_only,
    show_noncompliant,
  )

  const license_str = (license: string | undefined, url: string | undefined) =>
    url?.startsWith(`http`)
      ? `<a href="${url}" target="_blank" rel="noopener noreferrer" title="View license">${license}</a>`
      : `<span title="License file not available">${license}</span>`

  return (
    MODELS.filter(
      (model) => current_filter(model) && model.metrics?.discovery?.[discovery_set],
    )
      .map((model) => {
        const { license, metrics } = model
        const discover_metrics = metrics?.discovery?.[discovery_set]
        const is_compliant = model_is_compliant(model)

        // Get RMSD from geo_opt metrics if available, using the first symprec value
        const geo_opt_metrics = model.metrics?.geo_opt
        let rmsd = undefined
        if (geo_opt_metrics && typeof geo_opt_metrics === `object`) {
          // Try to find the first symprec key and get its RMSD
          const symprec_keys = Object.keys(geo_opt_metrics).filter((key) =>
            key.startsWith(`symprec=`),
          )
          if (symprec_keys.length > 0) {
            const symprec_key = symprec_keys[0]
            rmsd = get_geo_opt_property<number>(
              geo_opt_metrics,
              symprec_key.replace(`symprec=`, ``),
              `rmsd`,
            )
          }
        }

        // Get kappa from phonon metrics
        const phonons = metrics?.phonons
        const kappa =
          phonons && typeof phonons === `object` && `kappa_103` in phonons
            ? (phonons.kappa_103?.κ_SRME as number | undefined)
            : undefined

        const cps = calculate_cps(discover_metrics?.F1, rmsd, kappa, config)

        const targets = model.targets.replace(/_(.)/g, `<sub>$1</sub>`)
        const targets_str = `<span title="${targets_tooltips[model.targets]}">${targets}</span>`
        const row_style = show_noncompliant
          ? `border-left: 3px solid ${is_compliant ? compliant_clr : noncompliant_clr};`
          : null

        // Add model links
        const code_license = license?.code
          ? license_str(license.code, license.code_url)
          : `n/a`
        const checkpoint_license = license?.checkpoint
          ? license_str(license.checkpoint, license.checkpoint_url)
          : `n/a`

        const r_cut = model.hyperparams?.graph_construction_radius
        const r_cut_str = r_cut
          ? `<span data-sort-value="${r_cut}">${r_cut} Å</span>`
          : `n/a`

        return {
          Model: `<a title="Version: ${model.model_version}" href="/models/${model.model_key}" data-sort-value="${model.model_name}">${model.model_name}</a>`,
          CPS: cps,
          F1: discover_metrics?.F1,
          DAF: discover_metrics?.DAF,
          Prec: discover_metrics?.Precision,
          Acc: discover_metrics?.Accuracy,
          TPR: discover_metrics?.TPR,
          TNR: discover_metrics?.TNR,
          MAE: discover_metrics?.MAE,
          RMSE: discover_metrics?.RMSE,
          'R<sup>2</sup>': discover_metrics?.R2,
          'κ<sub>SRME</sub>': kappa,
          RMSD: rmsd,
          'Training Set': format_train_set(model.training_set, model),
          Params: `<span title="${pretty_num(model.model_params, `,`)}" trainable model parameters" data-sort-value="${model.model_params}">${pretty_num(model.model_params)}</span>`,
          Targets: targets_str,
          'Date Added': `<span title="${format_date(model.date_added)}" data-sort-value="${new Date(model.date_added).getTime()}">${model.date_added}</span>`,
          // Add Links as a special property
          Links: {
            paper: {
              url: model.paper || model.doi,
              title: `Read model paper`,
              icon: `<svg><use href="#icon-paper"></use></svg>`,
            },
            repo: {
              url: model.repo,
              title: `View source code`,
              icon: `<svg><use href="#icon-code"></use></svg>`,
            },
            pr_url: {
              url: model.pr_url,
              title: `View pull request`,
              icon: `<svg><use href="#icon-pull-request"></use></svg>`,
            },
            checkpoint: {
              url: model.checkpoint_url,
              title: `Download model checkpoint`,
              icon: `<svg><use href="#icon-download"></use></svg>`,
            },
            pred_files: { files: get_pred_file_urls(model), name: model.model_name },
          } as LinkData,
          'Checkpoint License': checkpoint_license,
          'Code License': code_license,
          'r<sub>cut</sub>': r_cut_str,
          row_style,
        }
      })
      // Sort by combined score (descending)
      .sort((row1, row2) => {
        const [score1, score2] = [row1[`CPS`], row2[`CPS`]]

        // Handle NaN values (they should be sorted to the bottom)
        const is_nan1 = score1 === null || isNaN(score1)
        const is_nan2 = score2 === null || isNaN(score2)
        if (is_nan1 && is_nan2) return 0
        if (is_nan1) return 1
        if (is_nan2) return -1

        return score2 - score1
      })
  )
}

// Calculate table cell background color based on its value and column config
export function calc_cell_color(
  val: number | null | undefined, // cell value
  all_values: (number | null | undefined)[], // all values in the column
  better: `higher` | `lower` | undefined, // sort direction
  color_scale: string | null = `interpolateViridis`, // color scale name
  scale_type: `linear` | `log` = `linear`, // scale type
): { bg: string | null; text: string | null } {
  // Skip color calculation for null values or if color_scale is null
  if (val === null || val === undefined || color_scale === null) {
    return { bg: null, text: null }
  }

  const numeric_vals = all_values.filter(
    (v): v is number => typeof v === `number` && v > 0, // Filter out non-positives for log scale
  )

  if (numeric_vals.length === 0) return { bg: null, text: null }

  const range = [min(numeric_vals) ?? 0, max(numeric_vals) ?? 1]

  // Reverse the range if lower values are better
  if (better === `lower`) range.reverse()

  // Use custom color scale if specified, otherwise fall back to viridis
  const scale_name = color_scale || `interpolateViridis`
  // Cast to ensure TypeScript recognizes it as a valid interpolator function
  const interpolator = (d3sc[scale_name as keyof typeof d3sc] ||
    d3sc.interpolateViridis) as (t: number) => string

  let color_scale_fn

  if (scale_type === `log` && range[0] > 0 && range[1] > 0) {
    // Use log scale for positive values
    color_scale_fn = scaleLog().domain(range).range([0, 1]).clamp(true)

    const normalized_val = color_scale_fn(val)
    const bg = interpolator(normalized_val)
    const text = choose_bw_for_contrast(null, bg)

    return { bg, text }
  } else {
    // Use sequential scale for linear mapping
    color_scale_fn = scaleSequential().domain(range).interpolator(interpolator)

    const bg = color_scale_fn(val)
    const text = choose_bw_for_contrast(null, bg)

    return { bg, text }
  }
}

export const CPS_COLUMN: HeatmapColumn = {
  label: `CPS`,
  tooltip: DEFAULT_CPS_CONFIG.description,
  style: `border-right: 1px solid black;`,
  format: `.3f`,
  better: `higher`,
}

export const METADATA_COLS: HeatmapColumn[] = [
  { label: `Model`, sticky: true, sortable: true, better: null },
  { label: `Training Set`, tooltip: `Size of and link to model training set` },
  { label: `Params`, tooltip: `Number of trainable model parameters` },
  { label: `Targets`, tooltip: `Target property used to train the model` },
  { label: `Date Added`, tooltip: `Submission date to the leaderboard` },
  {
    label: `Links`,
    tooltip: `Model resources: paper, code repository and submission pull request`,
    sortable: false,
  },
  {
    label: `r<sub>cut</sub>`,
    tooltip: `Graph construction radius in Ångströms (cutoff distance for creating edges in the graph)`,
    visible: false,
  },
  {
    label: `Checkpoint License`,
    tooltip: `Model checkpoint license`,
    visible: false,
  },
  {
    label: `Code License`,
    tooltip: `Model code license`,
    visible: false,
  },
]

export const DISCOVERY_METRICS: HeatmapColumn[] = [
  { label: `F1`, tooltip: `Harmonic mean of precision and recall` },
  { label: `DAF`, tooltip: `Discovery acceleration factor` },
  { label: `Prec`, tooltip: `Precision of classifying thermodynamic stability` },
  { label: `Acc`, tooltip: `Accuracy of classifying thermodynamic stability` },
  { label: `TPR`, tooltip: `True positive rate of classifying thermodynamic stability` },
  { label: `TNR`, tooltip: `True negative rate of classifying thermodynamic stability` },
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
]

export const PHONON_METRICS: HeatmapColumn[] = [
  {
    label: `κ<sub>SRME</sub>`,
    tooltip: `Symmetric relative mean error in predicted phonon mode contributions to thermal conductivity κ`,
    style: `border-left: 1px solid black;`,
  },
]

// Define geometry optimization metrics
export const GEO_OPT_METRICS: HeatmapColumn[] = [
  {
    label: `RMSD`,
    tooltip: `Root mean squared displacement between predicted and reference structures after relaxation`,
    style: `border-left: 1px solid black;`,
  },
  {
    label: `Energy Diff`,
    tooltip: `Mean absolute energy difference between predicted and reference structures`,
  },
  {
    label: `Force RMSE`,
    tooltip: `Root mean squared error of forces in predicted structures relative to reference`,
  },
  {
    label: `Stress RMSE`,
    tooltip: `Root mean squared error of stress in predicted structures relative to reference`,
  },
  {
    label: `Max Force`,
    tooltip: `Maximum force component in predicted structures after relaxation`,
  },
]

export const ALL_METRICS: HeatmapColumn[] = [
  CPS_COLUMN,
  ...DISCOVERY_METRICS,
  ...PHONON_METRICS,
  ...GEO_OPT_METRICS.slice(0, 1), // Only include RMSD by default, others can be toggled
]

// Define display labels and tooltips for metrics
export const metric_labels: Record<
  keyof ModelStats,
  { label: string; tooltip?: string; unit?: string }
> = {
  CPS: { label: `CPS`, tooltip: `Combined Performance Score` },
  Accuracy: { label: `Accuracy` },
  DAF: {
    label: `DAF`,
    tooltip: `Discovery Acceleration Factor`,
  },
  F1: { label: `F1 Score`, tooltip: `Harmonic mean of precision and recall` },
  MAE: { label: `MAE`, tooltip: `Mean Absolute Error`, unit: `eV / atom` },
  Precision: {
    label: `Precision`,
    tooltip: `Precision of classifying thermodynamic stability`,
  },
  R2: { label: `R<sup>2</sup>`, tooltip: `Coefficient of determination` },
  RMSE: { label: `RMSE`, tooltip: `Root Mean Squared Error`, unit: `eV / atom` },
  TNR: { label: `TNR`, tooltip: `True Negative Rate` },
  TPR: { label: `TPR`, tooltip: `True Positive Rate` },
  model_name: { label: `Model Name` },
  κ_SRME: {
    label: `κ<sub>SRME</sub>`,
    tooltip: `Symmetric relative mean error in predicted phonon mode contributions to thermal conductivity κ`,
  },
  // Add missing fields from ModelStats interface
  Recall: { label: `Recall`, tooltip: `Recall of classifying thermodynamic stability` },
  missing_preds: {
    label: `Missing Predictions`,
    tooltip: `Number of missing predictions`,
  },
  missing_percent: { label: `Missing %`, tooltip: `Percentage of missing predictions` },
  'Run Time (h)': { label: `Run Time`, tooltip: `Runtime in hours`, unit: `h` },
  GPUs: { label: `GPUs`, tooltip: `Number of GPUs used` },
  CPUs: { label: `CPUs`, tooltip: `Number of CPUs used` },
}

export const DISCOVERY_SET_LABELS: Record<
  DiscoverySet,
  { title: string; tooltip: string; link?: string }
> = {
  full_test_set: {
    title: `Full Test Set`,
    tooltip: `Metrics computed on the full test set including duplicate structure prototypes`,
  },
  unique_prototypes: {
    title: `Unique Prototypes`,
    tooltip: `Metrics computed only on ~215k unique structure prototypes in WBM determined by matching Aflow-style prototype strings.`,
    link: `https://github.com/janosh/matbench-discovery/blob/37baf7986f848/data/wbm/compile_wbm_test_set.py#L640-L654`,
  },
  most_stable_10k: {
    title: `10k Most Stable`,
    tooltip: `Metrics computed on the 10k structures predicted to be most stable (different for each model)`,
  },
}
