import DATASETS from '$data/datasets.yml'
import { format_date } from '$lib'
import {
  CDS_COMPONENTS,
  CDS_CONFIG,
  CMDS_CONFIG,
  CPS_CONFIG,
  type CdsPillar,
} from '$lib/combined-scores.svelte'
import { ACTIVE_MODELS, get_pred_file_urls } from '$lib/models.svelte'
import {
  ALL_METRICS,
  DISCOVERY_METRICS,
  DIATOMICS_METRICS,
  HYPERPARAMS,
  MD_METRICS,
  METADATA_COLS,
} from '$lib/labels'
import type { ModelMetadata, TargetType } from '$lib/schema/model'
import type { DiscoverySet, Label, ModelData } from '$lib/types'
import MODELINGS_TASKS from '$pkg/modeling-tasks.yml'
import { escape_html } from 'matterviz/utils'
import { format_num } from 'matterviz/labels'
import { is_invalid, type CellVal } from 'matterviz/table'

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

export const model_role_from_targets = (targets: TargetType) =>
  targets === `E`
    ? { label: `Energy predictor`, title: `Structure-to-energy predictor (no forces)` }
    : { label: `Interatomic potential`, title: `Force-capable interatomic potential` }

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

// Paths come from the fixed metric/metadata labels; reuse their parsed segments on redraw.
const data_paths = new Map<string, string[]>()
export function get_nested_value(model: ModelData, dotted_path: string): unknown {
  let keys = data_paths.get(dotted_path)
  if (!keys) {
    keys = dotted_path.split(`.`).filter(Boolean)
    data_paths.set(dotted_path, keys)
  }
  if (keys.length === 0) return undefined // empty path returns undefined, not the whole model
  let value: unknown = model

  for (const key of keys) {
    if (typeof value !== `object` || value === null) return undefined
    value = Reflect.get(value, key) // absent properties already yield undefined
  }

  return value
}

export function get_nested_number(
  model: ModelData,
  dotted_path: string,
): number | undefined {
  const value = get_nested_value(model, dotted_path)
  return typeof value === `number` ? value : undefined
}

export const is_finite_num = (value: unknown): value is number =>
  typeof value === `number` && Number.isFinite(value)

// Wrap a numeric value in a sortable span, or render `n/a` when undefined
const sortable_span = (value: number | undefined): string =>
  value === undefined ? `n/a` : `<span data-sort-value="${value}">${value}</span>`

// Build dot-separated data access path from a label, preferring `property` (actual
// data field name) over `key` when the two differ
export const label_data_path = (label: Label | undefined): string =>
  `${label?.path ?? ``}.${label?.property ?? label?.key ?? ``}`.replace(/^\./, ``)

export const metric_data_path = (
  label: Label,
  discovery_set: DiscoverySet = `unique_prototypes`,
): string =>
  label.path === DISCOVERY_METRICS.F1.path
    ? `metrics.discovery.${discovery_set}.${label.property ?? label.key}`
    : label_data_path(label)

export const metric_value = (
  model: ModelData,
  label: Label,
  discovery_set: DiscoverySet = `unique_prototypes`,
): unknown => get_nested_value(model, metric_data_path(label, discovery_set))

// Compose combined-score reasons without repeating the submission invitation.
export const missing_metric_reason = (model: ModelData, label: Label): string =>
  [...new Set(missing_metric_messages(model, label))]
    .toSorted(
      (left, right) =>
        Number(left.startsWith(`Contributions welcome`)) -
        Number(right.startsWith(`Contributions welcome`)),
    )
    .join(` `)

// Explain absent results using declared task status before inferring capability.
function missing_metric_messages(model: ModelData, label: Label): string[] {
  const is_cps = label.key === ALL_METRICS.CPS.key
  const raw_task =
    label.path?.split(`.`)[1] ??
    (label.key === MD_METRICS.md_time_multiplier.key
      ? `md`
      : label.key === DIATOMICS_METRICS.diatomics_time_multiplier.key
        ? `diatomics`
        : undefined)
  const task_key =
    raw_task && Object.hasOwn(MODELINGS_TASKS, raw_task)
      ? (raw_task as keyof NonNullable<ModelData['metrics']>)
      : undefined
  if (!task_key && !is_cps) {
    return [`${label.label.replaceAll(/<[^>]*>/g, ``)}: not reported.`]
  }

  const task_name = task_key ? MODELINGS_TASKS[task_key].label : `CPS`
  const task_data = task_key ? model.metrics?.[task_key] : undefined
  const detail = task_data?.reason ? ` ${task_data.reason}` : ``
  if (task_data?.status === `not_applicable`) {
    return [`${task_name}: unsupported.${detail}`]
  }
  const requires_forces =
    [`geo_opt`, `phonons`, `md`].includes(task_key ?? ``) ||
    (task_key === `diatomics` &&
      (label.key.includes(`force`) ||
        label.key === DIATOMICS_METRICS.diatomics_combined_score.key))
  if (model.targets === `E` && requires_forces) {
    return [`${task_name} requires forces; this model predicts only energies.${detail}`]
  }

  const invite = `Contributions welcome to add missing model predictions.`
  const property = label.property ?? label.key
  if (
    [`run_time_sec`, `max_rss_gb`, `max_gpu_mem_gb`].includes(property) ||
    label.key.endsWith(`time_multiplier`)
  ) {
    return [
      `${task_name}: ${label.key.endsWith(`time_multiplier`) ? `positive runtime` : property === `run_time_sec` ? `runtime` : `peak memory`} not reported.`,
      `Contributions welcome to add missing timing or memory data.`,
    ]
  }
  if (task_data?.status === `pending`) {
    return [`${task_name}: evaluation pending.${detail}`, invite]
  }
  if (task_data?.status === `not_available`) {
    return [`${task_name}: results unavailable.${detail}`, invite]
  }
  const has_results =
    task_data &&
    Object.keys(task_data).some((key) => key !== `status` && key !== `reason`)
  if (is_cps || (has_results && property === `combined_score`)) {
    const components: Label[] =
      is_cps || task_key === `md`
        ? Object.values(is_cps ? CPS_CONFIG : CMDS_CONFIG).filter(
            ({ weight }) => weight > 0,
          )
        : Object.entries(CDS_COMPONENTS).flatMap(([pillar, entries]) =>
            CDS_CONFIG[pillar as CdsPillar].weight > 0
              ? entries.map(({ key }) => ({
                  key,
                  label: key,
                  path: `metrics.diatomics`,
                  description: ``,
                }))
              : [],
          )
    const missing = components.filter((component) => {
      const value = metric_value(model, component)
      return (
        !is_finite_num(value) ||
        ((component.property ?? component.key) === `run_time_sec` && value <= 0)
      )
    })
    return missing.length
      ? missing.flatMap((component) => missing_metric_messages(model, component))
      : [`${label.label}: invalid components or weights.${detail}`]
  }
  return has_results
    ? [`${task_name}: incomplete results.${detail}`, invite]
    : [`${task_name}: not evaluated yet.${detail}`, invite]
}

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

type MetricsRowData = Pick<ModelData, `org_logos` | `authors`> & {
  model_key: string
  model_name: string
  Model: string
  CPS: ModelData[`CPS`]
  model: ModelData
  class?: string
  Links: Record<`paper` | `repo` | `pr_url` | `checkpoint`, string | null> & {
    pred_files: { files: { name: string; url: string }[]; name: string }
  }
  [key: string]: CellVal | ModelData[`org_logos` | `authors`]
}

const metadata_labels = [...Object.values(HYPERPARAMS), ...Object.values(METADATA_COLS)]

export function assemble_row_data(
  discovery_set: DiscoverySet,
  model_filter: (model: ModelData) => boolean,
  filter_matches: (model: ModelData) => boolean = () => true,
  models: ModelData[] = ACTIVE_MODELS, // injectable for tests
): MetricsRowData[] {
  const license_str = (license: string | undefined, url: string | null | undefined) =>
    url?.startsWith(`http`)
      ? `<a href="${url}" target="_blank" rel="noopener noreferrer" title="View license">${license}</a>`
      : `<span title="License file not available">${license}</span>`

  const filtered_models = models.filter(
    (model) => model_filter(model) && filter_matches(model),
  )

  const metric_num = (model: ModelData, label: Label) =>
    get_nested_number(model, label_data_path(label))
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
    const targets = model.targets.replaceAll(/_(?<char>.)/g, `<sub>$<char></sub>`)
    const targets_str = `<span title="${targets_tooltips[model.targets]}">${targets}</span>`

    const code_license = license?.code
      ? license_str(license.code, license.code_url)
      : `n/a`
    const checkpoint_license = license?.checkpoint
      ? license_str(license.checkpoint, license.checkpoint_url)
      : `n/a`

    const r_cut = model.hyperparams?.architecture?.graph_construction_radius
    const r_cut_str = r_cut ? `<span data-sort-value="${r_cut}">${r_cut} Å</span>` : `n/a`

    const { ase_optimizer, max_steps, max_force, cell_filter } =
      model.hyperparams?.evaluation ?? {}
    const { n_layers } = model.hyperparams?.architecture ?? {}
    const cell_filter_display =
      cell_filter && typeof cell_filter === `string`
        ? cell_filter.replace(/CellFilter$/, ``)
        : null
    const diatomics_metrics = metrics?.diatomics ?? null
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

    const row: MetricsRowData = {
      model_key: model.model_key,
      model_name: model.model_name,
      model,
      Model: `<a title="Version: ${model.model_version ?? `unknown`}" href="/models/${model.model_key}" data-sort-value="${model.model_name}">${model.model_name}</a>${model_exclusion_marker}`,
      ...Object.fromEntries(
        Object.values(ALL_METRICS).map((label) => [
          label.key,
          metric_value(model, label, discovery_set),
        ]),
      ),
      CPS: model.CPS,
      // computed after the spreads so they override the (pathless) spread entries
      [MD_METRICS.md_time_multiplier.key]: md_time_multiplier(model),
      [DIATOMICS_METRICS.diatomics_time_multiplier.key]: diatomics_time_multiplier(model),
      'Training Set': format_train_set(model.training_sets, model),
      [HYPERPARAMS.model_params.key]:
        `<span title="${format_num(model.model_params, `,`)} trainable model parameters" data-sort-value="${model.model_params}">${format_num(model.model_params)}</span>`,
      [HYPERPARAMS.ase_optimizer.key]: ase_optimizer ?? `n/a`,
      [HYPERPARAMS.max_steps.key]: sortable_span(max_steps),
      [HYPERPARAMS.max_force.key]: sortable_span(max_force),
      [HYPERPARAMS.cell_filter.key]: cell_filter_display
        ? `<span data-sort-value="${cell_filter}">${cell_filter_display}</span>`
        : `n/a`,
      [HYPERPARAMS.n_layers.key]: sortable_span(n_layers),
      Targets: targets_str,
      [METADATA_COLS.benchmark_added.key]:
        `<span title="${model.dates.benchmark_added ? format_date(model.dates.benchmark_added) : `Unknown`}" data-sort-value="${new Date(model.dates.benchmark_added ?? ``).getTime()}">${model.dates.benchmark_added ?? `n/a`}</span>`,
      Links: {
        paper: model.paper ?? model.doi,
        repo: model.repo,
        pr_url: model.pr_url,
        checkpoint: model.checkpoint_url,
        pred_files: { files: get_pred_file_urls(model), name: model.model_name },
      },
      [METADATA_COLS.checkpoint_license.label]: checkpoint_license,
      [METADATA_COLS.code_license.label]: code_license,
      [HYPERPARAMS.graph_construction_radius.key]: r_cut_str,
      org_logos: model.org_logos,
      authors: model.authors,
    }
    for (const label of metadata_labels) {
      if (is_invalid(row[label.key]) || row[label.key] === `n/a`) {
        row[label.key] =
          `<span data-title="${escape_html(missing_metric_reason(model, label))}">n/a</span>`
      }
    }
    return row
  })

  // Sort by combined performance score (descending)
  return all_metrics.toSorted((row1, row2) => {
    const score1 = row1.CPS ?? Number.NaN
    const score2 = row2.CPS ?? Number.NaN
    // NaN scores sort to the bottom
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

    const val_1 = get_nested_value(model_1, sort_by)
    const val_2 = get_nested_value(model_2, sort_by)

    // null/undefined/NaN sort last
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
