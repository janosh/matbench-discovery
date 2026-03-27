import { default as DATASETS } from '$data/datasets.yml'
import type { Author, ModelData } from '$lib/types'
import MODELINGS_TASKS from '$pkg/modeling-tasks.yml'
import { calculate_cps, CPS_CONFIG, type CpsConfig } from './combined_perf_score.svelte'
import { get_org_logo } from './labels'

export const MODEL_METADATA_PATHS = import.meta.glob<ModelData>(
  `$root/models/[^_]**/[^_]*.yml`,
  { eager: true, import: 'default' },
)

// Visually distinct color palette
export const MODEL_COLORS = [
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
      console.warn(`Training set ${data_name} not found in DATASETS`)
      continue
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
    .map(([key, metadata], index) => {
      // Assign color to each model for consistent coloring across plots
      const model_color = MODEL_COLORS[index % MODEL_COLORS.length]

      const sizes = calculate_training_sizes(metadata.training_set)

      // Get top affiliations with logos
      const affiliation_counts: Record<string, number> = {}
      const affiliation_data: Record<
        string,
        { name: string; id?: string; src?: string }
      > = {}

      for (const author of metadata.authors ?? ([] as Author[])) {
        if (!author.affiliation) continue

        const org_logo = get_org_logo(author.affiliation)
        const logo_key = org_logo?.id ?? org_logo?.src
        if (logo_key && org_logo) {
          affiliation_counts[logo_key] = (affiliation_counts[logo_key] || 0) + 1
          if (!(logo_key in affiliation_data)) {
            affiliation_data[logo_key] = org_logo
          }
        } else if (!import.meta.env.PROD) {
          // Only warn about missing logos in dev mode
          console.warn(`No logo found for affiliation: ${author.affiliation}`)
        }
      }

      const frequent_logos = Object.entries(affiliation_counts)
        .toSorted(([_key_a, count_a], [_key_b, count_b]) => count_b - count_a)
        .slice(0, 3)
        .map(([key]) => affiliation_data[key])

      return Object.assign({}, metadata, {
        dirname: key.split(`/`)[2],
        metadata_file: key.replace(/^..\//, ``),
        color: model_color,
        CPS: Number.NaN, // Initial CPS placeholder
        n_training_materials: sizes.total_materials,
        n_training_structures: sizes.total_structures,
        org_logos: frequent_logos,
      }) as ModelData
    }),
)

// Update CPSs of models based on current CPS weights
export function update_models_cps(models: ModelData[], cps_config: CpsConfig) {
  models.forEach((model: ModelData) => {
    // Extract required metrics for CPS calculation
    const discovery = model.metrics?.discovery
    const f1 =
      typeof discovery === `object` ? discovery?.[`unique_prototypes`]?.F1 : undefined
    const rmsd =
      model.metrics?.geo_opt && typeof model.metrics.geo_opt !== `string`
        ? model.metrics.geo_opt[`symprec=1e-5`]?.rmsd
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

// Compute compliant training sets from datasets.yml (datasets with compliant: true)
export const COMPLIANT_TRAINING_SETS: string[] = Object.entries(DATASETS)
  .filter(([_, val]) => typeof val === `object` && !Array.isArray(val) && val.compliant)
  .map(([key]) => key)

export function model_is_compliant(model: ModelData): boolean {
  if ((model.openness ?? `OSOD`) !== `OSOD`) return false

  return model.training_set.every((set) => COMPLIANT_TRAINING_SETS.includes(set))
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
