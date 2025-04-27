import { DATASETS, format_date, MODELS } from '$lib'
import type { TargetType } from '$lib/model-schema'
import { get_pred_file_urls, model_is_compliant } from '$lib/models.svelte'
import modeling_tasks from '$pkg/modeling-tasks.yml'
import { max, min } from 'd3-array'
import { scaleLog, scaleSequential } from 'd3-scale'
import * as d3sc from 'd3-scale-chromatic'
import { choose_bw_for_contrast, pretty_num } from 'elementari/labels'
import { GEO_OPT_SYMMETRY_METRICS, HYPERPARAMS, INFO_COLS, METRICS } from './labels'
import type { DiscoverySet, LinkData, ModelData } from './types'

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

// access (possibly deeply) nested metrics and parameters from ModelData objects
export function get_nested_value(
  model: ModelData, // model data to extract value from
  dotted_path: string, // dotted path to nested value to extract
): unknown {
  const keys = dotted_path.split(`.`).filter(Boolean) // remove empty parts
  let value: unknown = model

  for (const key of keys) {
    // Check if value is an object and has the key
    if (typeof value === `object` && value && key in value) {
      value = (value as Record<string, unknown>)[key]
    } else return undefined // Can't go deeper/property doesn't exist
  }

  return value
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

// combined filter function that respects both prediction type and compliance filters
export function make_combined_filter(
  model_filter: (model: ModelData) => boolean, // custom model filter function
  show_energy: boolean, // whether to show energy-only models
  show_compliant: boolean, // whether to show compliant models
  show_non_compliant: boolean, // whether to show non-compliant models
): (model: ModelData) => boolean {
  return (model: ModelData) => {
    if (!model_filter(model)) return false // Apply user-provided model_filter first

    // Filter energy-only models if not shown
    const is_energy_only = model.targets === `E`
    if (is_energy_only && !show_energy) return false

    // Filter noncompliant models if not shown
    const is_compliant = model_is_compliant(model)
    if (is_compliant && !show_compliant) return false
    if (!is_compliant && !show_non_compliant) return false

    return true
  }
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

// Calculate table data for the metrics table with combined scores
export function assemble_row_data(
  discovery_set: DiscoverySet, // discovery set to use for metrics
  model_filter: (model: ModelData) => boolean, // filter function for models
  show_energy_only: boolean, // show energy-only models
  show_non_compliant: boolean, // show non-compliant models
  show_compliant: boolean, // show compliant models
) {
  const current_filter = make_combined_filter(
    model_filter,
    show_energy_only,
    show_compliant,
    show_non_compliant,
  )

  const license_str = (license: string | undefined, url: string | undefined) =>
    url?.startsWith(`http`)
      ? `<a href="${url}" target="_blank" rel="noopener noreferrer" title="View license">${license}</a>`
      : `<span title="License file not available">${license}</span>`

  const filtered_models = MODELS.filter(
    (model) => current_filter(model) && model.metrics?.discovery?.[discovery_set],
  )

  const all_metrics = filtered_models.map((model) => {
    const { license, metrics } = model
    const discovery_metrics = metrics?.discovery?.[discovery_set]
    const is_compliant = model_is_compliant(model)
    const { RMSD, CPS } = METRICS

    // Get kappa from phonon metrics
    const phonons = metrics?.phonons
    const kappa =
      phonons && typeof phonons === `object` && `kappa_103` in phonons
        ? (phonons.kappa_103?.κ_SRME as number | undefined)
        : undefined

    const targets = model.targets.replace(/_(.)/g, `<sub>$1</sub>`)
    const targets_str = `<span title="${targets_tooltips[model.targets]}">${targets}</span>`
    const row_style = `border-left: 3px solid var(--${is_compliant ? `` : `non-`}compliant-color);`

    // Add model links
    const code_license = license?.code
      ? license_str(license.code, license.code_url)
      : `n/a`
    const checkpoint_license = license?.checkpoint
      ? license_str(license.checkpoint, license.checkpoint_url)
      : `n/a`

    const r_cut = model.hyperparams?.graph_construction_radius
    const r_cut_str = r_cut ? `<span data-sort-value="${r_cut}">${r_cut} Å</span>` : `n/a`

    return {
      Model: `<a title="Version: ${model.model_version}" href="/models/${model.model_key}" data-sort-value="${model.model_name}">${model.model_name}</a>`,
      CPS: model[CPS.key] as number | undefined,
      F1: discovery_metrics?.F1,
      DAF: discovery_metrics?.DAF,
      Prec: discovery_metrics?.Precision,
      Acc: discovery_metrics?.Accuracy,
      TPR: discovery_metrics?.TPR,
      TNR: discovery_metrics?.TNR,
      MAE: discovery_metrics?.MAE,
      RMSE: discovery_metrics?.RMSE,
      'R<sup>2</sup>': discovery_metrics?.R2,
      'κ<sub>SRME</sub>': kappa,
      RMSD: get_nested_value(model, `${RMSD.path}.${RMSD.key}`) as number | undefined,
      'Training Set': format_train_set(model.training_set, model),
      [HYPERPARAMS.model_params.short as string]:
        `<span title="${pretty_num(model.model_params, `,`)}" trainable model parameters" data-sort-value="${model.model_params}">${pretty_num(model.model_params)}</span>`,
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
      [INFO_COLS.checkpoint_license.label]: checkpoint_license,
      [INFO_COLS.code_license.label]: code_license,
      [HYPERPARAMS.graph_construction_radius.short as string]: r_cut_str,
      row_style,
      org_logos: model.org_logos,
      ...Object.fromEntries(
        Object.values(GEO_OPT_SYMMETRY_METRICS).map((col) => [
          col.label,
          get_nested_value(model, `${col.path}.${col.key}`) as number | undefined,
        ]),
      ),
    }
  })

  // Sort by combined performance score (descending)
  return all_metrics.sort((row1, row2) => {
    const [score1, score2] = [row1[`CPS`], row2[`CPS`]]

    // Handle undefined or null values (they should be sorted to the bottom)
    const is_nan1 = score1 == null || isNaN(score1)
    const is_nan2 = score2 == null || isNaN(score2)
    if (is_nan1 && is_nan2) return 0
    if (is_nan1) return 1
    if (is_nan2) return -1

    return score2 - score1
  })
}

// Sort models by a given metric (as dotted path sort_by) and order
export const sort_models =
  (sort_by: string, order: `asc` | `desc`) =>
  (model_1: ModelData, model_2: ModelData): number => {
    const sort_factor = order === `asc` ? -1 : 1

    // Special case for model_name sorting
    if (sort_by === `model_name`) {
      // For model_name, directly use localeCompare with sort_factor
      return sort_factor * model_1.model_name.localeCompare(model_2.model_name)
    }

    // Get values using the helper function for other metrics
    const val_1 = get_nested_value(model_1, sort_by)
    const val_2 = get_nested_value(model_2, sort_by)

    // Handle null, undefined, or NaN values by sorting last
    if (val_1 == null && val_2 == null) return 0
    if (val_1 == null || Number.isNaN(val_1)) return 1 // Always sort nulls/NaN to the end
    if (val_2 == null || Number.isNaN(val_2)) return -1 // Always sort nulls/NaN to the end

    if (typeof val_1 == `string` && typeof val_2 == `string`) {
      return sort_factor * (val_1 as string).localeCompare(val_2 as string)
    } else if (typeof val_1 == `number` && typeof val_2 == `number`) {
      // interpret runt_time==0 as infinity
      if (sort_by == `Run Time (h)`) {
        if (val_1 == 0) return -sort_factor
        if (val_2 == 0) return sort_factor
      }
      return sort_factor * (val_2 - val_1)
    } else {
      throw `Unexpected type '${val_1}' encountered sorting by key '${sort_by}'`
    }
  }
