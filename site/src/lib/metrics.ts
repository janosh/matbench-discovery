import { DATASETS, format_date, MODELS } from '$lib'
import {
  ALL_METRICS,
  DIATOMICS_METRICS,
  GEO_OPT_SYMMETRY_METRICS,
  HYPERPARAMS,
  MD_METRICS,
  METADATA_COLS,
} from '$lib/labels'
import type { ModelMetadata, ModelType, TargetType } from '$lib/schema/model'
import { get_pred_file_urls } from '$lib/models.svelte'
import type { DiscoverySet, Label, LinkData, ModelData } from '$lib/types'
import MODELINGS_TASKS from '$pkg/modeling-tasks.yml'
import { escape_html, format_num } from 'matterviz'

// Model target type descriptions
export const targets_tooltips: Record<TargetType, string> = {
  E: `Energy`,
  EF_G: `Energy with gradient-based forces`,
  EF_D: `Energy with direct forces`,
  EFS_G: `Energy with gradient-based forces and stress`,
  EFSH_G: `Energy with gradient-based forces, stress, and Hessian`,
  EFS_D: `Energy with direct forces and stress`,
  EFS_GM: `Energy with gradient-based forces, stress, and magmoms`,
  EFS_DM: `Energy with direct forces, stress, and magmoms`,
} as const

export const model_type_tooltips: Record<ModelType, string> = {
  GNN: `Graph Neural Network`,
  UIP: `Universal Interatomic Potential`,
  'BO-GNN': `Bayesian Optimization with Graph Neural Network`,
  Fingerprint: `Handcrafted feature-based model`,
  Transformer: `Attention-based transformer architecture`,
  RF: `Random Forest`,
} as const

export const openness_tooltips: Record<ModelMetadata[`openness`], string> = {
  OSOD: `Open source, open data`,
  OSCD: `Open source, closed data`,
  CSOD: `Closed source, open data`,
  CSCD: `Closed source, closed data`,
} as const

export const discovery_task_tooltips: Record<
  ModelMetadata[`train_task`] | ModelMetadata[`test_task`],
  string
> = {
  RP2RE: `relaxed prototype to relaxed energy`,
  RS2RE: `relaxed structure to relaxed energy`,
  S2E: `structure to energy`,
  S2RE: `structure to relaxed energy`,
  S2EF: `structure to energy, force`,
  S2EFS: `structure to energy, force, stress`,
  S2EFSM: `structure to energy, force, stress, magmoms`,
  IP2E: `initial prototype to energy`,
  IS2E: `initial structure to energy`,
  IS2RE: `initial structure to relaxed energy`,
  'IS2RE-SR': `initial structure to relaxed energy with structure relaxation`,
} as const

// Access (possibly deeply) nested metrics and parameters from ModelData objects
export function get_nested_value(model: ModelData, dotted_path: string): unknown {
  const keys = dotted_path.split(`.`).filter(Boolean) // Remove empty parts
  if (keys.length === 0) return undefined // Empty path returns undefined, not the whole model
  let value: unknown = model

  for (const key of keys) {
    // Check if value is an object and has the key
    if (typeof value === `object` && value && key in value) {
      value = Reflect.get(value, key) // dynamic lookup without weakening the surrounding type
    } else return undefined // Can't go deeper/property doesn't exist
  }

  return value
}

// Type-safe wrapper around get_nested_value that returns number | undefined
export function get_nested_number(
  model: ModelData,
  dotted_path: string,
): number | undefined {
  const value = get_nested_value(model, dotted_path)
  return typeof value === `number` ? value : undefined
}

// Type guard for finite numbers (excludes NaN, Infinity, and non-number values)
export const is_finite_num = (value: unknown): value is number =>
  typeof value === `number` && Number.isFinite(value)

// Wrap a numeric value in a sortable span, or render `n/a` when undefined
const sortable_span = (value: number | undefined): string =>
  value === undefined ? `n/a` : `<span data-sort-value="${value}">${value}</span>`

// Build dot-separated data access path from a label, preferring `property` (actual
// data field name) over `key` when the two differ
export const label_data_path = (label: Label | undefined): string =>
  `${label?.path ?? ``}.${label?.property ?? label?.key ?? ``}`.replace(/^\./, ``)

// Append "(higher|lower)=better" hint to a column tooltip where applicable
export function append_better_hint(col: Label, better = col.better): string {
  const description = col.description ?? ``
  if (better !== `higher` && better !== `lower`) return description
  return description ? `${description} (${better}=better)` : `${better}=better`
}

const all_higher_better_metrics = new Set(
  Object.values(MODELINGS_TASKS).flatMap(
    (model_task) => model_task.metrics.higher_is_better,
  ),
)

const all_lower_better_metrics = new Set(
  Object.values(MODELINGS_TASKS).flatMap(
    (model_task) => model_task.metrics.lower_is_better,
  ),
)

export function metric_better_as(metric: string): `higher` | `lower` | null {
  if (all_higher_better_metrics.has(metric)) return `higher`
  return all_lower_better_metrics.has(metric) ? `lower` : null
}

// Format training set information for display in the metrics table
export function format_train_set(model_train_sets: string[], model: ModelData): string {
  const { n_training_structures = 0, n_training_materials = 0 } = model

  const data_urls: Record<string, string> = {}
  const tooltip: string[] = []

  for (const data_name of model_train_sets) {
    if (!(data_name in DATASETS)) {
      throw new Error(`Training set ${data_name} not found in DATASETS`)
    }
    const { name, slug, n_structures, n_materials = n_structures } = DATASETS[data_name]
    data_urls[data_name] = `/data/${slug}`

    const structures_note =
      n_materials !== n_structures ? ` (${format_num(n_structures, `,`)} structures)` : ``
    tooltip.push(`${name}: ${format_num(n_materials, `,`)} materials${structures_note}`)
  }

  // render `_x` dataset-key suffixes as subscripts, e.g. ω_q -> ω<sub>q</sub>
  const sub = (key: string) =>
    key.replaceAll(/_(?<subscript>\w+)/g, `<sub>$<subscript></sub>`)
  const dataset_links = Object.entries(data_urls)
    .map(([key, href]) => `<a href="${href}">${sub(key)}</a>`)
    .join(`+`)
  const new_line = `&#013;` // Line break that works in title attribute
  const dataset_tooltip =
    tooltip.length > 1 ? `${new_line}• ${tooltip.join(`${new_line}• `)}` : ``

  const same_count = n_training_materials === n_training_structures
  const title = same_count
    ? `${format_num(n_training_materials, `,`)} materials in training set${new_line}${dataset_tooltip}`
    : `${format_num(n_training_materials, `,`)} materials in training set ` +
      `(${format_num(n_training_structures, `,`)} structures counting all DFT relaxation ` +
      `frames per material)${dataset_tooltip}`
  const sort_value = same_count
    ? n_training_materials
    : n_training_materials || n_training_structures
  const structure_count = same_count
    ? ``
    : ` <small>(${format_num(n_training_structures)})</small>`

  return (
    `<span title="${title}" data-sort-value="${sort_value}">` +
    `${format_num(n_training_materials)}${structure_count} ` +
    `<small>${dataset_links}</small></span>`
  )
}

// NB: cell background/text colors are computed by matterviz's HeatmapTable internally
// (calc_cell_color in matterviz/table) — no local color logic needed

// Calculate table data for the metrics table with combined scores
export function assemble_row_data(
  discovery_set: DiscoverySet,
  model_filter: (model: ModelData) => boolean,
  filter_matches: (model: ModelData) => boolean = () => true,
  models: ModelData[] = MODELS, // injectable for tests
) {
  const license_str = (license: string | undefined, url: string | undefined) =>
    url?.startsWith(`http`)
      ? `<a href="${url}" target="_blank" rel="noopener noreferrer" title="View license">${license}</a>`
      : `<span title="License file not available">${license}</span>`

  const filtered_models = models.filter(
    (model) => model_filter(model) && filter_matches(model),
  )

  const { RMSD } = ALL_METRICS
  // label_data_path prefers label.property over label.key, so columns whose row key
  // must differ from the YAML field (e.g. the two run_time_sec columns) resolve too
  const metric_num = (model: ModelData, label: Label) =>
    get_nested_number(model, label_data_path(label))
  const metric_columns = <Labels extends Record<keyof Labels, Label>>(
    model: ModelData,
    labels: Labels,
  ) =>
    Object.fromEntries(
      Object.values<Label>(labels).map((label) => [label.key, metric_num(model, label)]),
    )
  const finite_positive = (value: unknown): value is number =>
    is_finite_num(value) && value > 0
  // Slowdown columns: wall time relative to the fastest model in the current
  // filtered view (roster-dependent, so computed here rather than stored on models)
  const time_multiplier = (run_time_label: Label) => {
    const fastest = Math.min(
      ...filtered_models
        .map((model) => metric_num(model, run_time_label))
        .filter(finite_positive),
    )
    return (model: ModelData) => {
      const run_time = metric_num(model, run_time_label)
      return finite_positive(run_time) ? run_time / fastest : undefined
    }
  }
  const md_time_multiplier = time_multiplier(MD_METRICS.md_run_time_sec)
  const diatomics_time_multiplier = time_multiplier(
    DIATOMICS_METRICS.diatomics_run_time_sec,
  )
  const all_metrics = filtered_models.map((model) => {
    const { license, metrics } = model
    const discovery_metrics =
      typeof metrics?.discovery === `object`
        ? metrics.discovery[discovery_set]
        : undefined
    const targets = model.targets.replaceAll(/_(?<char>.)/g, `<sub>$<char></sub>`)
    const targets_str = `<span title="${targets_tooltips[model.targets]}">${targets}</span>`

    // Add model links
    const code_license = license?.code
      ? license_str(license.code, license.code_url)
      : `n/a`
    const checkpoint_license = license?.checkpoint
      ? license_str(license.checkpoint, license.checkpoint_url)
      : `n/a`

    const r_cut = model.hyperparams?.graph_construction_radius
    const r_cut_str = r_cut ? `<span data-sort-value="${r_cut}">${r_cut} Å</span>` : `n/a`

    // Get geometry optimization hyperparameters
    const { ase_optimizer, max_steps, max_force, cell_filter, n_layers } =
      model.hyperparams ?? {}
    const cell_filter_display =
      cell_filter && typeof cell_filter === `string`
        ? cell_filter.replace(/CellFilter$/, ``)
        : null
    const diatomics_metrics =
      typeof metrics?.diatomics === `object` ? metrics.diatomics : null
    const excluded_formula_reasons = diatomics_metrics?.excluded_formula_reasons ?? {}
    // group excluded formulas by reason for a compact tooltip like
    // "Diatomics metrics exclude A-A, B-B due to <reason>; C-C due to <other>"
    let model_exclusion_marker = ``
    if (Object.keys(excluded_formula_reasons).length > 0) {
      // manual grouping instead of Map.groupBy, which is newer than Vite's default
      // browser baseline and would throw at runtime in e.g. Safari < 17.4
      const formulas_by_reason = new Map<string, string[]>()
      for (const [formula, reason] of Object.entries(excluded_formula_reasons)) {
        const group = formulas_by_reason.get(reason) ?? []
        group.push(formula)
        formulas_by_reason.set(reason, group)
      }
      const exclusion_note = escape_html(
        `Diatomics metrics exclude ${[...formulas_by_reason]
          .map(
            ([reason, formulas]) =>
              `${formulas.join(`, `)}${reason ? ` due to ${reason}` : ``}`,
          )
          .join(`; `)}`,
      )
      model_exclusion_marker =
        `<span title="${exclusion_note}" aria-label="${exclusion_note}">` +
        `<span aria-hidden="true">*</span></span>`
    }

    return {
      model_name: model.model_name,
      Model: `<a title="Version: ${model.model_version}" href="/models/${model.model_key}" data-sort-value="${model.model_name}">${model.model_name}</a>${model_exclusion_marker}`,
      CPS: model.CPS,
      F1: discovery_metrics?.F1,
      DAF: discovery_metrics?.DAF,
      Precision: discovery_metrics?.Precision,
      Recall: discovery_metrics?.Recall,
      Accuracy: discovery_metrics?.Accuracy,
      TPR: discovery_metrics?.TPR,
      TNR: discovery_metrics?.TNR,
      MAE: discovery_metrics?.MAE,
      RMSE: discovery_metrics?.RMSE,
      R2: discovery_metrics?.R2,
      [ALL_METRICS.κ_SRME.key]: metric_num(model, ALL_METRICS.κ_SRME),
      [ALL_METRICS.κ_SRE.key]: metric_num(model, ALL_METRICS.κ_SRE),
      [RMSD.key]: metric_num(model, RMSD),
      ...metric_columns(model, MD_METRICS),
      ...metric_columns(model, DIATOMICS_METRICS),
      // computed after the spreads so they override the (pathless) spread entries
      [MD_METRICS.md_time_multiplier.key]: md_time_multiplier(model),
      [DIATOMICS_METRICS.diatomics_time_multiplier.key]: diatomics_time_multiplier(model),
      'Training Set': format_train_set(model.training_set, model),
      [HYPERPARAMS.model_params.key]:
        `<span title="${format_num(model.model_params, `,`)} trainable model parameters" data-sort-value="${model.model_params}">${format_num(model.model_params)}</span>`,
      [HYPERPARAMS.ase_optimizer.key]: ase_optimizer ?? `n/a`,
      [HYPERPARAMS.max_steps.key]: sortable_span(max_steps),
      [HYPERPARAMS.max_force.key]: sortable_span(max_force),
      [HYPERPARAMS.cell_filter.key]: cell_filter_display
        ? `<span data-sort-value="${cell_filter}">${cell_filter_display}</span>`
        : `n/a`,
      [HYPERPARAMS.n_layers.key]: sortable_span(n_layers),
      ...metric_columns(model, GEO_OPT_SYMMETRY_METRICS),
      Targets: targets_str,
      [METADATA_COLS.date_added.key]:
        `<span title="${format_date(model.date_added)}" data-sort-value="${new Date(model.date_added).getTime()}">${model.date_added}</span>`,
      Links: {
        paper: {
          url: model.paper || model.doi,
          title: `Read model paper`,
          icon: `Paper`,
        },
        repo: { url: model.repo, title: `View source code`, icon: `Code` },
        pr_url: { url: model.pr_url, title: `View pull request`, icon: `PullRequest` },
        checkpoint: {
          url: model.checkpoint_url,
          title: `Download model checkpoint`,
          icon: `Download`,
        },
        pred_files: { files: get_pred_file_urls(model), name: model.model_name },
      } satisfies LinkData,
      [METADATA_COLS.checkpoint_license.label]: checkpoint_license,
      [METADATA_COLS.code_license.label]: code_license,
      [HYPERPARAMS.graph_construction_radius.key]: r_cut_str,
      org_logos: model.org_logos,
      authors: model.authors,
    }
  })

  // Sort by combined performance score (descending)
  return all_metrics.toSorted((row1, row2) => {
    const score1 = row1.CPS ?? Number.NaN
    const score2 = row2.CPS ?? Number.NaN
    // Handle missing or NaN values (they should be sorted to the bottom)
    if (Number.isNaN(score1)) return Number.isNaN(score2) ? 0 : 1
    return Number.isNaN(score2) ? -1 : score2 - score1
  })
}

// Sort models by a given metric (as dotted path sort_by) and order
export const sort_models =
  (sort_by: string, order: `asc` | `desc`) =>
  (model_1: ModelData, model_2: ModelData): number => {
    const sort_factor = order === `asc` ? 1 : -1

    // Special case for Model sorting (by model_name): asc = alphabetical A->Z
    if (sort_by === `Model`) {
      return sort_factor * model_1.model_name.localeCompare(model_2.model_name)
    }

    // Get values using the helper function for other metrics
    const val_1 = get_nested_value(model_1, sort_by)
    const val_2 = get_nested_value(model_2, sort_by)

    // Handle null, undefined, or NaN values by sorting last
    const sorts_last = (val: unknown) =>
      val == null || (typeof val === `number` && Number.isNaN(val))
    if (sorts_last(val_1) && sorts_last(val_2)) return 0
    if (sorts_last(val_1)) return 1
    if (sorts_last(val_2)) return -1

    if (typeof val_1 === `string` && typeof val_2 === `string`) {
      return sort_factor * val_1.localeCompare(val_2)
    }
    if (typeof val_1 === `number` && typeof val_2 === `number`) {
      // Interpret run_time === 0 as infinity
      if (sort_by === `Run Time`) {
        if (val_1 === 0 && val_2 === 0) return 0
        if (val_1 === 0) return sort_factor
        if (val_2 === 0) return -sort_factor
      }
      return sort_factor * (val_1 - val_2)
    }
    throw new TypeError(
      `Unexpected type '${typeof val_1}' encountered sorting by key '${sort_by}'`,
    )
  }
