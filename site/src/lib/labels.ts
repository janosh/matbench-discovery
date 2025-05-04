import MODELINGS_TASKS from '$pkg/modeling-tasks.yml'
import type { DiscoverySet, Metric } from './types'

export const RMSD_BASELINE = 0.15 // baseline for poor performance given worst performing model at time of writing is M3GNet at 0.1117

// Helper function to format scientific notation with superscript
// used e.g. for symprec in geo_opt metrics
export const format_power_ten = (text: string): string => {
  return text
    .replace(/(\d+(?:\.\d+)?)e\+?(-?\d+)/gi, (_, base, exponent) => {
      return `${base}×10<sup>${exponent}</sup>`
    })
    .replace(`1×10`, `10`)
}

export const sub_sup_to_tspan = (text: string): string => {
  return text
    .replaceAll(`<sup>`, `<tspan baseline-shift='0.4em' font-size='0.8em'>`)
    .replaceAll(`</sup>`, `</tspan>`)
}

export const DISCOVERY_METRICS: Record<string, Metric> = {
  Accuracy: {
    key: `Accuracy`,
    short: `Acc`,
    label: `Accuracy`,
    description: `Accuracy of classifying thermodynamic stability`,
    better: `higher`,
    path: `metrics.discovery.unique_prototypes`,
  },
  F1: {
    key: `F1`,
    short: `F1`,
    label: `F1 Score`,
    path: `metrics.discovery.unique_prototypes`,
    description: `Harmonic mean of precision and recall for stable/unstable material classification (discovery task)`,
    range: [0, 1],
    better: `higher`,
  },
  DAF: {
    key: `DAF`,
    label: `DAF`,
    description: `Discovery Acceleration Factor measuring how much better ML models classify thermodynamic stability compared to random guessing`,
    better: `higher`,
    path: `metrics.discovery.unique_prototypes`,
  },
  Precision: {
    key: `Precision`,
    short: `Prec`,
    label: `Precision`,
    description: `Precision of classifying thermodynamic stability`,
    path: `metrics.discovery.unique_prototypes`,
    better: `higher`,
  },
  Recall: {
    key: `Recall`,
    label: `Recall`,
    description: `Recall of classifying thermodynamic stability`,
    path: `metrics.discovery.unique_prototypes`,
    better: `higher`,
    visible: false,
  },
  TNR: {
    key: `TNR`,
    label: `TNR`,
    description: `True Negative Rate`,
    path: `metrics.discovery.unique_prototypes`,
    better: `higher`,
  },
  TPR: {
    key: `TPR`,
    label: `TPR`,
    description: `True Positive Rate`,
    path: `metrics.discovery.unique_prototypes`,
    better: `higher`,
  },
  MAE: {
    key: `MAE`,
    label: `MAE`,
    description: `Mean Absolute Error`,
    path: `metrics.discovery.unique_prototypes`,
    better: `lower`,
  },
  R2: {
    key: `R2`,
    label: `R<sup>2</sup>`,
    description: `Coefficient of determination`,
    path: `metrics.discovery.unique_prototypes`,
    better: `higher`,
  },
  RMSE: {
    key: `RMSE`,
    label: `RMSE`,
    description: `Root Mean Squared Error`,
    path: `metrics.discovery.unique_prototypes`,
    better: `lower`,
  },
} as const

export const METADATA_COLS: Record<string, Metric> = {
  model_name: {
    label: `Model`,
    key: `model_name`,
    description: `Model name`,
    sticky: true,
    sortable: true,
    better: null,
  },
  training_set: {
    label: `Training Set`,
    key: `training_set`,
    description: `Size of and link to model training set`,
  },
  targets: {
    label: `Targets`,
    key: `targets`,
    description: `Target property used to train the model`,
  },
  date_added: {
    label: `Date Added`,
    key: `date_added`,
    format: `%b %y`,
    description: `Submission date to the leaderboard`,
  },
  links: {
    label: `Links`,
    key: `links`,
    description: `Model resources: paper, code repository and submission pull request`,
    sortable: false,
  },
  r_cut: {
    label: `r<sub>cut</sub>`,
    key: `r_cut`,
    description: `Graph construction radius in Ångströms (cutoff distance for creating edges in the graph)`,
    visible: false,
  },
  n_training_materials: {
    label: `Number of Training Materials`,
    key: `n_training_materials`,
    description: `Number of training materials`,
    format: `~s`,
  },
  n_training_structures: {
    label: `Number of Training Structures`,
    key: `n_training_structures`,
    description: `Number of training structures`,
    format: `~s`,
  },
  checkpoint_license: {
    label: `Checkpoint License`,
    key: `checkpoint_license`,
    description: `Model checkpoint license`,
    visible: false,
  },
  code_license: {
    label: `Code License`,
    key: `code_license`,
    description: `Model code license`,
    visible: false,
  },
  missing_preds: {
    key: `missing_preds`,
    label: `Missing Predictions`,
    description: `Number of missing predictions`,
    visible: false,
  },
  missing_percent: {
    key: `missing_percent`,
    label: `Missing %`,
    description: `Percentage of missing predictions`,
    visible: false,
  },
  'Run Time (h)': {
    key: `run_time_h`,
    label: `Run Time`,
    description: `Runtime in hours`,
    visible: false,
    better: `lower`,
  },
  org: {
    key: `org`,
    label: `Org`,
    sortable: false,
    description: `Most common author affiliations`,
    style: `text-align: center; max-width: 2em; transform: scale(1.2);`,
    visible: true,
    better: null,
  },
} as const

export const HYPERPARAMS: Record<string, Metric> = {
  model_params: {
    label: `Number of model parameters`,
    key: `model_params`,
    short: `Params`,
    description: `Number of trainable model parameters`,
    format: `~s`,
  },
  graph_construction_radius: {
    label: `Graph construction radius r<sub>cut</sub>`,
    key: `graph_construction_radius`,
    short: `r<sub>cut</sub>`,
    path: `hyperparams`,
    description: `Graph construction radius in Ångströms (cutoff distance for creating edges in the graph)`,
  },
  max_force: {
    label: `Max force`,
    key: `max_force`,
    path: `hyperparams`,
    description: `Maximum force in eV/Å`,
  },
  max_steps: {
    label: `Max relaxation steps`,
    key: `max_steps`,
    path: `hyperparams`,
    description: `Maximum number of steps`,
  },
  batch_size: {
    label: `Batch size`,
    key: `batch_size`,
    path: `hyperparams`,
    description: `Batch size`,
  },
  epochs: {
    label: `Training epochs`,
    key: `epochs`,
    path: `hyperparams`,
    description: `Number of training epochs`,
  },
  n_layers: {
    label: `Number of layers`,
    key: `n_layers`,
    path: `hyperparams`,
    description: `Number of (usually message passing) layers`,
  },
  learning_rate: {
    label: `Learning rate`,
    key: `learning_rate`,
    path: `hyperparams`,
    description: `Learning rate`,
  },
  max_neighbors: {
    label: `Max number of neighbors during graph construction`,
    key: `max_neighbors`,
    path: `hyperparams`,
    description: `Maximum number of neighbors`,
  },
  n_estimators: {
    label: `Number of estimators`,
    key: `n_estimators`,
    path: `hyperparams`,
    description: `Number of estimators`,
  },
} as const

export type MetricKey = keyof typeof ALL_METRICS

export const DATASET_METADATA_COLS: Record<string, Metric> = {
  name: { label: `Name`, key: `name`, sticky: true },
  structures: {
    key: `n_structures`,
    label: `Number of Structures`,
    short: `Structures`,
    description: `Number of structures in the dataset. Any system with atomic positions and energy/force/stress labels is counted as a structure incl. successive ionic steps in MD/geometry optimization trajectories.`,
    better: `higher`,
    scale_type: `log`,
    format: `.3s`,
  },
  materials: {
    key: `n_materials`,
    label: `Number of Materials`,
    short: `Materials`,
    better: `higher`,
    scale_type: `log`,
    format: `.3s`,
    description: `Number of unique materials/prototypes in the dataset.`,
  },
  created: {
    key: `created`,
    label: `Created`,
    description: `Date the dataset was created/started`,
  },
  open: { label: `Open`, key: `Open`, style: `text-align: center;` },
  static: {
    label: `Static`,
    key: `Static`,
    style: `text-align: center;`,
    description: `Whether the dataset is static (fixed version) or dynamic (continuously updated).`,
  },
  license: { key: `license`, label: `License` },
  method: { key: `method`, label: `Method`, style: `max-width: 5em;` },
  api: {
    label: `API`,
    key: `API`,
    description: `API docs (OPTIMADE or native)`,
    sortable: false,
  },
  links: { key: `links`, label: `Links`, sortable: false },
} as const

export const GEO_OPT_SYMMETRY_METRICS: Record<string, Metric> = Object.fromEntries(
  [`1e-2`, `1e-5`]
    .flatMap(
      (symprec) =>
        [
          [`symmetry_match`, `Σ<sub>=</sub>`, `higher`, symprec],
          [`symmetry_decrease`, `Σ<sub>↓</sub>`, `lower`, symprec],
          [`symmetry_increase`, `Σ<sub>↑</sub>`, undefined, symprec],
        ] as const,
    )
    .map(([key, label, better, symprec]) => [
      `${key}_${symprec}`,
      {
        key,
        symprec,
        path: `metrics.geo_opt.symprec=${symprec}`,
        short: `${label} ${format_power_ten(symprec)}`,
        label: `${label} (symprec=${format_power_ten(symprec)})`,
        svg_label: `${label} (symprec=${sub_sup_to_tspan(format_power_ten(symprec))})`
          .replace(`Σ<sub>`, `Σ<tspan baseline-shift='-0.4em' font-size='0.8em'>`)
          .replace(`</sub>`, `</tspan>`),
        description: `Fraction of structures where ML and DFT ground state have matching spacegroup at ${format_power_ten(symprec)} symprec`,
        better,
        format: `~%`,
        visible: false,
      },
    ]),
)

export const ALL_METRICS: Record<string, Metric> = {
  // Dynamic metrics
  CPS: {
    key: `CPS`,
    short: `CPS`,
    label: `Combined Performance Score`,
    description: `Combined Performance Score averages discovery (F1), structure optimization (RMSD), and phonon performance (κ<sub>SRME</sub>) according to user-defined weights`,
    range: [0, 1],
    better: `higher`,
    format: `.3f`,
  },
  ...DISCOVERY_METRICS,
  // Phonon metrics
  κ_SRME: {
    key: `κ_SRME`,
    label: `κ<sub>SRME</sub>`,
    svg_label: `κ<tspan baseline-shift='-0.4em' font-size='0.8em'>SRME</tspan>`,
    description: `Symmetric relative mean error in predicted phonon mode contributions to thermal conductivity κ`,
    path: `metrics.phonons.kappa_103`,
    better: `lower`,
  },
  // Geometry optimization metrics
  RMSD: {
    key: `rmsd`,
    path: `metrics.geo_opt.symprec=1e-2`,
    label: `RMSD`,
    range: [0, RMSD_BASELINE],
    better: `lower`,
    description: `Root mean squared displacement between predicted and reference structures after relaxation`,
    style: `border-left: 1px solid black;`,
  },
  ...DISCOVERY_METRICS,
  ...GEO_OPT_SYMMETRY_METRICS,
} as const

export const DISCOVERY_SET_LABELS: Record<
  DiscoverySet,
  { label: string; description: string; link?: string }
> = {
  full_test_set: {
    label: `Full Test Set`,
    description: `Metrics computed on the full 257k WBM test set including duplicate structure prototypes`,
  },
  unique_prototypes: {
    label: `Unique Prototypes`,
    description: `Metrics computed only on ~215k unique structure prototypes in WBM determined by matching Aflow-style prototype strings.`,
    link: `https://github.com/janosh/matbench-discovery/blob/37baf7986f848/data/wbm/compile_wbm_test_set.py#L640-L654`,
  },
  most_stable_10k: {
    label: `10k Most Stable`,
    description: `Metrics computed on the 10k structures predicted to be most stable (different for each model)`,
  },
} as const

export const PROPERTY_LABELS: Record<string, string> = Object.fromEntries(
  Object.values({ ...ALL_METRICS, ...METADATA_COLS, ...HYPERPARAMS }).map((prop) => [
    prop.key,
    prop.label,
  ]),
)

// Formats a property path for display in UI components
export function format_property_path(path: string): string {
  // Split path into components
  let parts = path
    .split(`.`)
    .filter((part) => ![`metrics`, `kappa_103`].includes(part) && part)

  // remove symprec value preceding rmsd
  if (parts.at(-1)?.toUpperCase() === `RMSD`) {
    parts = parts.filter((part) => !part.includes(`symprec`))
  }

  // Default formatting for other dotted paths
  return parts
    .map((part) => {
      const pretty_label = CATEGORY_LABELS[part] ?? PROPERTY_LABELS[part]
      if (pretty_label) return pretty_label
      return format_power_ten(part).replaceAll(`_`, ` `)
    })
    .join(` > `)
}

export const CATEGORY_LABELS = Object.fromEntries(
  Object.entries({ ...MODELINGS_TASKS, ...DISCOVERY_SET_LABELS, ...HYPERPARAMS }).map(
    ([key, task]) => [key, task.label],
  ),
)

// TODO maybe remove get_format() since unused
// Determines appropriate string format for displaying a set of numerical values
// based on their characteristics (magnitude, precision, etc.).
export function get_format(values: number[]): string {
  if (!values.length) return `.1f`

  const avg = values.reduce((sum, val) => sum + val, 0) / values.length
  const max = Math.max(...values)
  const min = Math.min(...values)

  // Check if values are in plausible date timestamps after Jan 1, 2000 (946684800000) and before Jan 1, 2050 (2524608000000)
  const vals_are_dates = min > 946684800000 && max < 2524608000000
  if (vals_are_dates) return `%b %y`

  // Format selection logic based on data characteristics
  if (max > 10000 || avg > 1000) return `.1s`
  if (Math.abs(avg) > 0 && Math.abs(avg) < 0.01) return `.5f`
  if (max - min > 1000) return `.2s`
  if (values.every((val) => Math.abs(val - Math.round(val)) < 1e-6)) return `d`

  return `.2f`
}

export const to_title = (str: string) => str.charAt(0).toUpperCase() + str.slice(1)
export const title_case = (str: string) =>
  str.replaceAll(`_`, ` `).split(` `).map(to_title).join(` `)

// Map of author affiliations in model YAMLs to SVG icons (either inline symbol ID
// or external file path under /static/logos/) and full affiliation names for tooltips. Each item can have SVG ID from app.html OR src path under /static/logos/.
export const org_logos = {
  'AI for Science Institute, Beijing': `/logos/beijing-ai-for-science-institute.svg`,
  'Argonne National Laboratory': `/logos/argonne-national-lab.svg`,
  'Chinese Academy of Sciences': `/logos/chinese-academy-of-sciences.svg`,
  'Cornell University': `/logos/cornell-university.svg`,
  'Deep Principle': `/logos/deep-principle.svg`,
  DeePMD: `/logos/deepmd.svg`,
  'FAIR at Meta': `icon-logo-meta`,
  'Google DeepMind': `/logos/deepmind.svg`,
  'ICAMS, Ruhr University Bochum': `/logos/icams-bochum.svg`,
  'Incheon National University': `/logos/incheon-national-university.svg`,
  'Massachusetts Institute of Technology': `/logos/mit.svg`,
  'Microsoft Research': `icon-logo-microsoft`,
  'National Institute of Standards and Technology': `/logos/nist.svg`,
  'Northwestern University': `/logos/northwestern-university.svg`,
  'Orbital Materials': `/logos/orbital-materials.svg`,
  'Seoul National University': `/logos/seoul-national-university.svg`,
  'Tsinghua University': `/logos/tsinghua-university.svg`,
  'UC San Diego': `/logos/uc-san-diego.svg`,
  'UC Berkeley': `/logos/uc-berkeley.svg`,
  'University of Cambridge': `/logos/cambridge-university.svg`,
  'University of Florida': `/logos/university-of-florida.svg`,
  'University of Minnesota': `/logos/university-of-minnesota.svg`,
  'University of Texas at Austin': `/logos/university-of-texas-austin.svg`,
  'Beijing Information Science and Technology University': `/logos/beijing-information-science-and-technology-university.svg`,
} as const

// Attempts to find a matching logo data (ID or src) and name for a given affiliation string.
// Performs a case-insensitive search for keywords defined in org_logos.
// Returns object with logo name and either SVG ID or src path, or undefined if no match.
export function get_org_logo(
  affiliation: string,
): { name: string; id?: string; src?: string } | undefined {
  if (!affiliation) return undefined

  for (const [key, logo] of Object.entries(org_logos)) {
    // Check if lowercased affiliation string includes lowercased org key
    if (affiliation.includes(key)) {
      if (logo.startsWith(`/logos/`)) return { name: key, src: logo }
      else return { name: key, id: logo }
    }
  }
}
