import { metric_labels } from './metrics'

// Determines the appropriate string format for displaying a set of numerical values
// based on their characteristics (magnitude, precision, etc.).
export function get_format(values: number[]): string {
  if (!values.length) return `.1f`

  const avg = values.reduce((sum, val) => sum + val, 0) / values.length
  const max = Math.max(...values)
  const min = Math.min(...values)

  // Check if values are in plausible date timestamps after Jan 1, 2000 (946684800000) and before Jan 1, 2050 (2524608000000)
  const vals_are_dates = min > 946_684_800_000 && max < 2_524_608_000_000
  if (vals_are_dates) return `%b %y`

  // Format selection logic based on data characteristics
  if (max > 10000 || avg > 1000) return `.1s`
  if (Math.abs(avg) > 0 && Math.abs(avg) < 0.01) return `.5f`
  if (max - min > 1000) return `.2s`
  if (values.every((val) => Math.abs(val - Math.round(val)) < 1e-6)) return `d`

  return `.2f`
}

// Category mapping for property paths
export const CATEGORY_LABELS: Record<string, string> = {
  discovery: `Discovery`,
  phonons: `Phonons`,
  geo_opt: `Geometry Optimization`,
  hyperparams: `Hyperparams`,
}

// Special property mapping for direct properties
export const PROPERTY_LABELS: Record<string, string> = {
  model_params: `Model Parameters`,
  n_estimators: `Number of Estimators`,
  date_added: `Date Added`,
  n_training_materials: `Number of Training Materials`,
  n_training_structures: `Number of Training Structures`,
  graph_construction_radius: `Graph Construction Radius r<sub>cut</sub>`,
  max_neighbors: `Max Neighbors`,
  max_force: `Max Force (eV/Å)`,
  max_steps: `Max Relaxation Steps`,
  learning_rate: `Learning Rate`,
  batch_size: `Batch Size`,
  epochs: `Training Epochs`,
  n_layers: `Number of Layers`,
  // Add metric names for clearer labels
  rmsd: `RMSD`,
  κ_SRME: `κ<sub>SRME</sub>`,
  CPS: `Combined Performance Score`,
}

// Helper function to format scientific notation with superscript
export const format_scientific_notation = (text: string): string => {
  return text
    .replace(/(\d+(?:\.\d+)?)e-(\d+)/gi, (_, base, exponent) => {
      return `${base}×10<sup>-${exponent}</sup>`
    })
    .replace(/(\d+(?:\.\d+)?)e\+?(\d+)/gi, (_, base, exponent) => {
      return `${base}×10<sup>${exponent}</sup>`
    })
    .replace(`1×10`, `10`)
}

// Formats a property path for display in UI components
export function format_property_path(path: string): string {
  // Handle direct properties without dots
  if (!path.includes(`.`)) return PROPERTY_LABELS[path] || path

  // Split path into components
  const parts = path.split(`.`)

  // Special case for category-based paths
  if (parts.length >= 2) {
    const [category, ...rest] = parts
    const category_label = CATEGORY_LABELS[category] || category

    // Handle discovery metrics (discovery.set.metric)
    if (category === `discovery` && parts.length === 3) {
      const [_, set, metric] = parts
      // Format discovery set name
      const formatted_set = set
        .split(`_`)
        .map((word) => word.charAt(0).toUpperCase() + word.slice(1))
        .join(` `)

      // Get metric label from metric_labels if available
      const metric_label =
        metric_labels[metric as keyof typeof metric_labels]?.label || metric

      return `${category_label} > ${formatted_set} > ${metric_label}`
    }

    // Handle hyperparameters (hyperparams.param)
    if (category === `hyperparams` && parts.length === 2) {
      const param = parts[1]
      return `${category_label} > ${
        PROPERTY_LABELS[param] ||
        param
          .split(`_`)
          .map((word) => word.charAt(0).toUpperCase() + word.slice(1))
          .join(` `)
      }`
    }

    // Handle geo_opt metrics with symprec pattern
    if (category === `geo_opt` && parts.length === 3 && parts[1].startsWith(`symprec=`)) {
      const [_category, _symprec, metric] = parts
      if (metric === `rmsd`) return `${category_label} > RMSD`
    }

    // Handle phonons and other categories
    if ([`phonons`].includes(category)) {
      if (category === `phonons` && parts.length === 3 && parts[1] === `kappa_103`) {
        return `${category_label} > κ<sub>SRME</sub>`
      }
      return `${category_label} > ${rest.join(` > `)}`
    }

    // Default formatting for other dotted paths
    return parts
      .map((part, idx) => {
        if (idx === 0) return CATEGORY_LABELS[part] || part
        if (idx === parts.length - 1) {
          return (
            metric_labels[part as keyof typeof metric_labels]?.label || part.toUpperCase()
          ) // For metrics like RMSD
        }
        // Apply scientific notation formatting to any part that might contain it
        return format_scientific_notation(part)
      })
      .join(` > `)
  }

  // Fallback for any other format
  return path
}
