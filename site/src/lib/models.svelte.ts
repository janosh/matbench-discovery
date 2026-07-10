import { default as DATASETS } from '$data/datasets.yml'
import type { DiscoverySet, ModelData } from '$lib/types'
import MODELINGS_TASKS from '$pkg/modeling-tasks.yml'
import type { CpsConfig } from '$lib/combined-scores.svelte'
import {
  calculate_cps,
  CPS_CONFIG,
  CDS_CONFIG,
  update_models_cds,
  CMDS_CONFIG,
  update_models_cmds,
} from './combined-scores.svelte'
import { get_org_logo, type OrgLogo } from './labels'
import { UrlTableFilters } from './url-state.svelte'

export const MODEL_METADATA_PATHS = import.meta.glob<ModelData>(
  `$root/models/[^_]**/[^_]*.yml`,
  { eager: true, import: 'default' },
)

// Visually distinct color palette
const MODEL_COLORS = [
  `#4285F4`, // Blue
  `#EA4335`, // Red
  `#FBBC05`, // Yellow
  `#34A853`, // Green
  `#8A2BE2`, // Blueviolet
  `#FF7F50`, // Coral
  `#1E90FF`, // Dodgerblue
  `#FF1493`, // Deeppink
  `#32CD32`, // Limegreen
  `#FF8C00`, // Darkorange
  `#9370DB`, // Mediumpurple
  `#3CB371`, // Mediumseagreen
  `#DC143C`, // Crimson
  `#6495ED`, // Cornflowerblue
  `#FFD700`, // Gold
  `#8B008B`, // Darkmagenta
  `#00CED1`, // Darkturquoise
  `#FF4500`, // Orangered
  `#2E8B57`, // Seagreen
  `#BA55D3`, // Mediumorchid
]

// Calculate the total number of materials and structures in a model's training set
export function calculate_training_sizes(model_train_sets: string[] = []): {
  total_materials: number
  total_structures: number
} {
  let total_materials = 0
  let total_structures = 0

  for (const data_name of model_train_sets) {
    if (!(data_name in DATASETS)) {
      throw new Error(`Training set ${data_name} not found in DATASETS`)
    }
    const { n_structures, n_materials = n_structures } = DATASETS[data_name]
    total_materials += n_materials
    total_structures += n_structures
  }

  return { total_materials, total_structures }
}

export const MODELS = $state(
  Object.entries(MODEL_METADATA_PATHS)
    .filter(
      // Ignore models with status != completed (the default status)
      ([_key, metadata]) => (metadata?.status ?? `complete`) === `complete`,
    )
    .map(([key, metadata], index): ModelData => {
      // Assign color to each model for consistent coloring across plots
      const model_color = MODEL_COLORS[index % MODEL_COLORS.length]

      const sizes = calculate_training_sizes(metadata.training_set)

      const org_logos: OrgLogo[] = []
      for (const author of metadata.authors ?? []) {
        const org_logo = author.affiliation ? get_org_logo(author.affiliation) : undefined

        if (author.affiliation && !org_logo && !import.meta.env.PROD) {
          // Only warn about missing logos in dev mode
          console.warn(`No logo found for affiliation: ${author.affiliation}`)
        }

        if (org_logo && !org_logos.some((logo) => logo.name === org_logo.name)) {
          org_logos.push(org_logo)
        }
      }

      return {
        ...metadata,
        dirname: key.split(`/`)[2],
        metadata_file: key.replace(/^..\//, ``),
        color: model_color,
        CPS: Number.NaN, // Initial CPS placeholder
        n_training_materials: sizes.total_materials,
        n_training_structures: sizes.total_structures,
        org_logos,
      }
    }),
)

// Update CPSs of models based on current CPS weights
export function update_models_cps(models: ModelData[], cps_config: CpsConfig) {
  models.forEach((model: ModelData) => {
    // Extract required metrics for CPS calculation
    const discovery = model.metrics?.discovery
    const f1 =
      typeof discovery === `object` ? discovery?.unique_prototypes?.F1 : undefined
    // use symprec=1e-2 to match the RMSD column path declared in ALL_METRICS.RMSD,
    // so the displayed RMSD is the same value that feeds into CPS
    const rmsd =
      model.metrics?.geo_opt && typeof model.metrics.geo_opt !== `string`
        ? model.metrics.geo_opt[`symprec=1e-2`]?.rmsd
        : undefined
    const kappa =
      model.metrics?.phonons && typeof model.metrics.phonons !== `string`
        ? model.metrics.phonons.kappa_103?.κ_SRME !== undefined
          ? Number(model.metrics.phonons.kappa_103.κ_SRME)
          : undefined
        : undefined

    // Calculate and update CPS
    model.CPS = calculate_cps(f1, rmsd, kappa, cps_config) ?? Number.NaN
  })
}

// Calculate initial CPS for all models
update_models_cps(MODELS, CPS_CONFIG)

// Calculate initial CMDS (combined MD score) for all models. Computed on the fly
// (never stored in model YAMLs) so it tracks the current formula and live reweighting.
update_models_cmds(MODELS, CMDS_CONFIG)

// Calculate initial CDS (combined diatomics score) for all models, same on-the-fly
// semantics as CPS/CMDS
update_models_cds(MODELS, CDS_CONFIG)

// All dataset keys used by at least one model's training_set, in datasets.yml
// declaration order — the roster for the table's training-data filter dropdown
export const ALL_TRAINING_SETS: string[] = Object.keys(DATASETS).filter((key) =>
  MODELS.some((model) => model.training_set.some((dataset) => dataset === key)),
)

// table filter (training data + openness + targets + heatmap) with the dataset roster
export const make_table_filters = (): UrlTableFilters =>
  new UrlTableFilters(ALL_TRAINING_SETS)

// Find the model with the highest F1 score on the given discovery set among the
// passed models (pre-filter to the table cohort so the landing summary describes
// the same models the table shows). Returns null if none qualifies.
export function find_best_model(
  models: ModelData[],
  discovery_set: DiscoverySet = `full_test_set`,
): ModelData | null {
  let best: ModelData | null = null
  let best_f1 = -Infinity

  for (const model of models) {
    const discovery = model.metrics?.discovery
    const f1 = typeof discovery === `object` ? discovery?.[discovery_set]?.F1 : undefined
    if (typeof f1 !== `number` || Number.isNaN(f1)) continue
    if (f1 > best_f1) [best, best_f1] = [model, f1]
  }
  return best
}

export function get_pred_file_urls(model: ModelData) {
  // Get all pred_file_url from model.metrics
  const files: { name: string; url: string }[] = []

  function find_pred_files(obj: object, parent_key = ``) {
    if (!obj || typeof obj !== `object`) return

    for (const [key, val] of Object.entries(obj)) {
      if (key === `pred_file_url` && val && typeof val === `string`) {
        // Look up the label by traversing the MODELINGS_TASKS hierarchy
        const pretty_label = get_label_for_key_path(parent_key)
        files.push({ name: pretty_label, url: val })
      } else if (typeof val === `object`) {
        find_pred_files(val, key)
      }
    }
  }

  // Recursively look up labels in the MODELINGS_TASKS object
  function get_label_for_key_path(key_path: string): string {
    const tasks = MODELINGS_TASKS
    if (key_path in tasks) return tasks[key_path].label

    // Check if it's a subtask by searching all tasks
    for (const task_value of Object.values(tasks)) {
      if (task_value.subtasks?.[key_path]) {
        return task_value.subtasks[key_path].label
      }
    }

    return key_path // Default to key itself if no label is found
  }

  if (model.metrics) find_pred_files(model.metrics)
  return files
}
