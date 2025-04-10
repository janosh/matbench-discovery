import { DATASETS, MODELS, get_pred_file_urls, model_is_compliant } from '$lib'
import type { TargetType } from '$lib/model-schema'
import { discovery as discovery_config } from '$pkg/modeling-tasks.yml'
import { max, min } from 'd3-array'
import { scaleLog, scaleSequential } from 'd3-scale'
import * as d3sc from 'd3-scale-chromatic'
import { choose_bw_for_contrast, pretty_num } from 'elementari/labels'
import type { DiscoveryMetricsSet } from './model-schema.d.ts'
import type {
  CombinedMetricConfig,
  DiscoverySet,
  HeatmapColumn,
  LinkData,
  ModelData,
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

  // Check for κ_SRME in phonons
  if (
    metric_key === `κ_SRME` &&
    model.metrics?.phonons &&
    typeof model.metrics.phonons !== `string`
  ) {
    const value =
      model.metrics.phonons.kappa_103?.κ_SRME ?? model.metrics.phonons[metric_key]
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

// Determines if a specific metric is considered "better" when lower
export function is_lower_better(metric_key: string): boolean {
  // First check discovery metrics
  if (discovery_config.metrics.lower_is_better.includes(metric_key)) {
    return true
  }

  // Check phonons metrics
  if (metric_key === `κ_SRME` || metric_key === `RMSD`) {
    return true
  }

  // For any other metric, assume higher is better
  return false
}

// Format date string into human-readable format
export function format_date(date: string, options?: Intl.DateTimeFormatOptions): string {
  return new Date(date).toLocaleDateString(undefined, {
    year: `numeric`,
    month: `short`,
    day: `numeric`,
    ...options,
  })
}

// Format training set information for display in the metrics table
export function format_train_set(model_train_sets: string[]): string {
  let [total_structs, total_materials] = [0, 0]
  const data_urls: Record<string, string> = {}
  const tooltip: string[] = []

  for (const train_set of model_train_sets) {
    if (!(train_set in DATASETS)) {
      console.warn(`Training set ${train_set} not found in DATASETS`)
      continue
    }
    const training_set_info = DATASETS[train_set]
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

  const dataset_links = data_links.join(`+`)
  const new_line = `&#013;` // line break that works in title attribute
  const dataset_tooltip =
    tooltip.length > 1 ? `${new_line}• ${tooltip.join(new_line + `• `)}` : ``

  let title = `${pretty_num(total_materials, `,`)} materials in training set${new_line}${dataset_tooltip}`
  let train_size_str = `<span title="${title}" data-sort-value="${total_materials}">${pretty_num(total_materials)} <small>${dataset_links}</small></span>`

  if (total_materials !== total_structs) {
    title =
      `${pretty_num(total_materials, `,`)} materials in training set ` +
      `(${pretty_num(total_structs, `,`)} structures counting all DFT relaxation ` +
      `frames per material)${dataset_tooltip}`

    train_size_str =
      `<span title="${title}" data-sort-value="${total_materials}">` +
      `${pretty_num(total_materials)} <small>(${pretty_num(total_structs)})</small> ` +
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
          'Training Set': format_train_set(model.training_set),
          Params: `<span title="${pretty_num(model.model_params, `,`)}" trainable model parameters" data-sort-value="${model.model_params}">${pretty_num(model.model_params)}</span>`,
          Targets: targets_str,
          'Date Added': `<span title="${format_date(model.date_added)}" data-sort-value="${new Date(model.date_added).getTime()}">${model.date_added}</span>`,
          // Add Links as a special property
          Links: {
            paper: {
              url: model.paper || model.doi,
              title: `Read model paper`,
              icon: `📄`,
            },
            repo: { url: model.repo, title: `View source code`, icon: `📦` },
            pr_url: { url: model.pr_url, title: `View pull request`, icon: `🔗` },
            checkpoint: {
              url: model.checkpoint_url,
              title: `Download model checkpoint`,
              icon: `💾`,
            },
            pred_files: { files: get_pred_file_urls(model), name: model.model_name },
          } as LinkData,
          'Checkpoint License': license_str(license?.checkpoint, license?.checkpoint_url),
          'Code License': license_str(license?.code, license?.code_url),
          row_style,
        }
      })
      // Sort by combined score (descending)
      .sort((row1, row2) => {
        const score1 = row1[`CPS`]
        const score2 = row2[`CPS`]

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
  const interpolator = d3sc[scale_name as keyof typeof d3sc] || d3sc.interpolateViridis

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

export const RMSD_BASELINE = 0.15 // baseline for poor performance given worst performing model at time of writing is M3GNet at 0.1117

export const DEFAULT_CPS_CONFIG: CombinedMetricConfig = {
  label: `CPS`,
  name: `Combined Performance Score`,
  key: `cps`,
  description: `Combined Performance Score averages discovery (F1), structure optimization (RMSD), and phonon performance (κ<sub>SRME</sub>) according to user-defined weights`,
  parts: {
    F1: {
      path: `discovery.unique_prototypes.F1`,
      label: `F1`,
      description: `F1 score for stable/unstable material classification (discovery task)`,
      weight: 0.5,
      range: [0, 1],
      better: `higher`,
    },
    kappa_SRME: {
      path: `phonons.kappa_103.κ_SRME`,
      label: `κ<sub>SRME</sub>`,
      svg_label: `κ<tspan baseline-shift='-0.4em' font-size='0.8em'>SRME</tspan>`,
      description: `Symmetric relative mean error of predicted lattice thermal conductivity`,
      weight: 0.4,
      range: [0, 2],
      better: `lower`,
    },
    RMSD: {
      path: `discovery.unique_prototypes.RMSD`,
      label: `RMSD`,
      description: `Root mean square displacement for crystal structure optimization`,
      weight: 0.1,
      range: [0, RMSD_BASELINE],
      better: `lower`,
    },
  },
} as const

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
    label: `Checkpoint License`,
    tooltip: `Model checkpoint license`,
    sortable: true,
    visible: false,
  },
  {
    label: `Code License`,
    tooltip: `Model code license`,
    sortable: true,
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

export const CPS_COLUMN: HeatmapColumn = {
  label: `CPS`,
  tooltip: DEFAULT_CPS_CONFIG.description,
  style: `border-right: 1px solid black;`,
  format: `.3f`,
  better: `higher`,
}

export const ALL_METRICS: HeatmapColumn[] = [
  CPS_COLUMN,
  ...DISCOVERY_METRICS,
  ...PHONON_METRICS,
  ...GEO_OPT_METRICS.slice(0, 1), // Only include RMSD by default, others can be toggled
]

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

// F1 score is between 0-1 where higher is better (no normalization needed)
function normalize_f1(value: number | undefined): number {
  if (value === undefined || isNaN(value)) return 0
  return value // Already in [0,1] range
}

// RMSD is lower=better, with current models in the range of ~0.01-0.25 Å
// We invert this so that better performance = higher score
function normalize_rmsd(value: number | undefined): number {
  if (value === undefined || isNaN(value)) return 0

  // Fixed reference points for RMSD (in Å)
  const excellent = 0 // Perfect performance (atoms in exact correct positions)

  // Linear interpolation between fixed points with clamping
  // Inverse mapping since lower RMSD is better
  if (value <= excellent) return 1.0
  if (value >= RMSD_BASELINE) return 0.0
  return (RMSD_BASELINE - value) / (RMSD_BASELINE - excellent)
}

// kappa_SRME is symmetric relative mean error, with range [0,2] by definition
// Lower values are better (0 is perfect)
function normalize_kappa_srme(value: number | undefined): number {
  if (value === undefined || isNaN(value)) return 0

  // Simple linear normalization from [0,2] to [1,0]
  // No clamping needed as SRME is bounded by definition
  return Math.max(0, 1 - value / 2)
}

// Calculate a combined score using normalized metrics weighted by importance factors.
// This uses fixed normalization reference points to ensure score stability when new models are added.

// Normalization reference points:
// - F1 score for discovery already in [0,1] range, higher is better
// - RMSD Root mean square displacement in range 0 (perfect) to RMSD_BASELINE, lower is better
// - κ_SRME symmetric relative mean error for lattice thermal conductivity,
//    range [0,2] linearly mapped to [1,0], lower is better
export function calculate_cps(
  f1: number | undefined,
  rmsd: number | undefined,
  kappa: number | undefined,
  config: CombinedMetricConfig, // weights for each metric
): number | null {
  // Find weights from config by metric names
  const { F1, RMSD, kappa_SRME } = config.parts

  // Check if any metrics with non-zero weights are missing
  if (
    (F1.weight > 0 && (f1 === undefined || isNaN(f1))) ||
    (RMSD.weight > 0 && (rmsd === undefined || isNaN(rmsd))) ||
    (kappa_SRME.weight > 0 && (kappa === undefined || isNaN(kappa)))
  ) {
    return null
  }

  // Skip the calculation if all weights are zero
  const total_weight = F1.weight + RMSD.weight + kappa_SRME.weight
  if (total_weight === 0) {
    return 0
  }

  // Calculate weighted sum
  let weighted_sum = 0

  // Add F1 contribution if available and weighted
  if (f1 !== undefined && !isNaN(f1) && F1.weight > 0) {
    weighted_sum += normalize_f1(f1) * F1.weight
  }

  // Add RMSD contribution if available and weighted
  if (rmsd !== undefined && !isNaN(rmsd) && RMSD.weight > 0) {
    weighted_sum += normalize_rmsd(rmsd) * RMSD.weight
  }

  // Add kappa contribution if available and weighted
  if (kappa !== undefined && !isNaN(kappa) && kappa_SRME.weight > 0) {
    weighted_sum += normalize_kappa_srme(kappa) * kappa_SRME.weight
  }

  // Return weighted average
  return weighted_sum / total_weight
}
