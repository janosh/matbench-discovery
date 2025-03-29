import TRAINING_SETS from '$data/training-sets.yml'
import { MODEL_METADATA, get_pred_file_urls, model_is_compliant } from '$lib'
import { pretty_num } from 'elementari'
import { calculate_combined_score } from './metrics'
import type { TargetType } from './model-schema.d.ts'
import type {
  CombinedMetricConfig,
  DiscoverySet,
  LinkData,
  ModelData,
  TableData,
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

// format date string into human-readable format
export function format_long_date(date: string): string {
  return new Date(date).toLocaleDateString(undefined, {
    weekday: `long`,
    year: `numeric`,
    month: `long`,
    day: `numeric`,
  })
}

// format training set information for display in the metrics table
export function format_train_set(model_training_sets: string[]): string {
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

// safely access nested geometry optimization properties
export function get_geo_opt_property<T>(
  geo_opt: unknown,
  symprec: string,
  property: string,
): T | undefined {
  if (!geo_opt || typeof geo_opt !== `object`) return undefined

  const symprec_key = `symprec=${symprec}`
  const metrics = geo_opt[symprec_key as keyof typeof geo_opt]

  if (!metrics || typeof metrics !== `object`) return undefined

  return metrics[property as keyof typeof metrics] as T
}

/**
 * Creates a combined filter function that respects both prediction type
 * and compliance filters.
 *
 * @param model_filter - The base model filter function
 * @param show_energy - Whether to show energy-only models
 * @param show_noncomp - Whether to show non-compliant models
 * @returns Combined filter function for models
 */
export function create_combined_filter(
  model_filter: (model: ModelData) => boolean,
  show_energy: boolean,
  show_noncomp: boolean,
): (model: ModelData) => boolean {
  return (model: ModelData) => {
    if (!model_filter(model)) {
      return false // Apply user-provided model_filter first
    }

    // Filter energy-only models if not shown
    const is_energy_only = model.targets === `E`
    if (is_energy_only && !show_energy) {
      return false
    }

    // Filter noncompliant models if not shown
    const is_compliant = model_is_compliant(model)
    if (!is_compliant && !show_noncomp) {
      return false
    }

    return true
  }
}

// calculate table data for the metrics table with combined scores
export function calculate_metrics_data(
  discovery_set: DiscoverySet,
  model_filter: (model: ModelData) => boolean,
  show_energy_only: boolean,
  show_noncompliant: boolean,
  config: CombinedMetricConfig,
  compliant_clr: string = `#4caf50`,
  noncompliant_clr: string = `#4682b4`,
): TableData {
  // Get the current filter with current state values
  const current_filter = create_combined_filter(
    model_filter,
    show_energy_only,
    show_noncompliant,
  )

  return (
    MODEL_METADATA.filter(
      (model) => current_filter(model) && model.metrics?.discovery?.[discovery_set],
    )
      .map((model) => {
        const metrics = model.metrics?.discovery?.[discovery_set]
        const is_compliant = model_is_compliant(model)

        // Get RMSD from geo_opt metrics if available, using the first symprec value
        const geo_opt_metrics = model.metrics?.geo_opt
        let rmsd = undefined
        if (geo_opt_metrics && typeof geo_opt_metrics === `object`) {
          // Try to find the first symprec key and get its RMSD
          const symprec_keys = Object.keys(geo_opt_metrics).filter((k) =>
            k.startsWith(`symprec=`),
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
        const phonons = model.metrics?.phonons
        const kappa =
          phonons && typeof phonons === `object` && `kappa_103` in phonons
            ? (phonons.kappa_103?.κ_SRME as number | undefined)
            : undefined

        // Calculate combined score
        const cps = calculate_combined_score(metrics?.F1, rmsd, kappa, config)

        const targets = model.targets.replace(/_(.)/g, `<sub>$1</sub>`)
        const targets_str = `<span title="${targets_tooltips[model.targets]}">${targets}</span>`
        const row_style = show_noncompliant
          ? `border-left: 3px solid ${is_compliant ? compliant_clr : noncompliant_clr};`
          : null

        return {
          Model: `<a title="Version: ${model.model_version}" href="/models/${model.model_key}">${model.model_name}</a>`,
          CPS: cps,
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
          RMSD: rmsd,
          'Training Set': format_train_set(model.training_set),
          Params: `<span title="${pretty_num(model.model_params, `,`)}" trainable model parameters" data-sort-value="${model.model_params}">${pretty_num(model.model_params)}</span>`,
          Targets: targets_str,
          'Date Added': `<span title="${format_long_date(model.date_added)}" data-sort-value="${new Date(model.date_added).getTime()}">${model.date_added}</span>`,
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
          row_style,
        }
      })
      // Sort by combined score (descending)
      .sort((row1, row2) => {
        // Handle NaN values (they should be sorted to the bottom)
        const score1 = row1[`CPS`]
        const score2 = row2[`CPS`]

        if (isNaN(score1) && isNaN(score2)) return 0
        if (isNaN(score1)) return 1
        if (isNaN(score2)) return -1

        return score2 - score1
      })
  )
}
