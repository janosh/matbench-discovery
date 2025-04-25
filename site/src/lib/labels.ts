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

// keys, display labels and descriptions for metrics
export const METRICS: Record<string, Metric> = {
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
  // Discovery metrics
  Accuracy: {
    key: `Accuracy`,
    short: `Acc`,
    label: `Accuracy`,
    description: `Accuracy of classifying thermodynamic stability`,
    better: `higher`,
    path: `metrics.discovery.unique_prototypes`,
  },
  DAF: {
    key: `DAF`,
    label: `DAF`,
    description: `Discovery Acceleration Factor measuring how much better ML models classify thermodynamic stability compared to random guessing`,
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
    short: `Rec`,
    label: `Recall`,
    description: `Recall of classifying thermodynamic stability`,
    path: `metrics.discovery.unique_prototypes`,
    better: `higher`,
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
} as const

export const METADATA_COLS: Record<string, Metric> = {
  model_name: {
    label: `Model`,
    key: `model_name`,
    sticky: true,
    sortable: true,
    better: null,
  },
  training_set: {
    label: `Training Set`,
    key: `training_set`,
    description: `Size of and link to model training set`,
  },
  model_params: {
    label: `Params`,
    key: `model_params`,
    format: `~s`,
    description: `Number of trainable model parameters`,
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
    label: `Training Materials`,
    key: `n_training_materials`,
    description: `Number of training materials`,
  },
  n_training_structures: {
    label: `Training Structures`,
    key: `n_training_structures`,
    description: `Number of training structures`,
  },
} as const

export const HYPERPARAMS: Record<string, Metric> = {
  graph_construction_radius: {
    label: `r<sub>cut</sub>`,
    key: `r_cut`,
    path: `hyperparams`,
    description: `Graph construction radius in Ångströms (cutoff distance for creating edges in the graph)`,
  },
  max_force: {
    label: `Max Force`,
    key: `max_force`,
    path: `hyperparams`,
    description: `Maximum force in eV/Å`,
  },
  max_steps: {
    label: `Max Steps`,
    key: `max_steps`,
    path: `hyperparams`,
    description: `Maximum number of steps`,
  },
  batch_size: {
    label: `Batch Size`,
    key: `batch_size`,
    path: `hyperparams`,
    description: `Batch size`,
  },
  epochs: {
    label: `Epochs`,
    key: `epochs`,
    path: `hyperparams`,
    description: `Number of training epochs`,
  },
  n_layers: {
    label: `Layers`,
    key: `n_layers`,
    path: `hyperparams`,
    description: `Number of (usually message passing) layers`,
  },
} as const

export const INFO_COLS: Record<string, Metric> = {
  // Metadata
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
  },
  missing_percent: {
    key: `missing_percent`,
    label: `Missing %`,
    description: `Percentage of missing predictions`,
  },
  'Run Time (h)': {
    key: `run_time_h`,
    label: `Run Time`,
    description: `Runtime in hours`,
  },
} as const

export type MetricKey = keyof typeof METRICS

export const DATASET_METADATA_COLS: Record<string, Metric> = {
  title: { key: `title`, label: `Title`, sticky: true },
  structures: {
    key: `n_structures`,
    label: `Number of Structures`,
    short: `Structures`,
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
  },
  created: {
    key: `created`,
    label: `Created`,
    description: `Date the dataset was created/started`,
  },
  open: { key: `open`, label: `Open` },
  license: { key: `license`, label: `License` },
  method: { key: `method`, label: `Method` },
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
        label: `${label} ${format_power_ten(symprec)}`,
        description: `Fraction of structures where ML and DFT ground state have matching spacegroup at ${format_power_ten(symprec)} symprec`,
        better,
        format: `.1%`,
        visible: false,
      },
    ]),
)

export const DISCOVERY_METRICS = Object.fromEntries(
  (
    [`F1`, `DAF`, `Precision`, `Accuracy`, `TPR`, `TNR`, `MAE`, `RMSE`, `R2`] as const
  ).map((key) => [key, METRICS[key]]),
)

export const ALL_METRICS: Record<string, Metric> = {
  [METRICS.CPS.key]: METRICS.CPS,
  ...DISCOVERY_METRICS,
  // Phonon metrics
  [METRICS.κ_SRME.key]: METRICS.κ_SRME,
  // Geometry optimization metrics
  [METRICS.RMSD.key]: METRICS.RMSD,
  ...GEO_OPT_SYMMETRY_METRICS,
} as const

export const DISCOVERY_SET_LABELS: Record<
  DiscoverySet,
  { title: string; description: string; link?: string }
> = {
  full_test_set: {
    title: `Full Test Set`,
    description: `Metrics computed on the full test set including duplicate structure prototypes`,
  },
  unique_prototypes: {
    title: `Unique Prototypes`,
    description: `Metrics computed only on ~215k unique structure prototypes in WBM determined by matching Aflow-style prototype strings.`,
    link: `https://github.com/janosh/matbench-discovery/blob/37baf7986f848/data/wbm/compile_wbm_test_set.py#L640-L654`,
  },
  most_stable_10k: {
    title: `10k Most Stable`,
    description: `Metrics computed on the 10k structures predicted to be most stable (different for each model)`,
  },
} as const

// Formats a property path for display in UI components
export function format_property_path(path: string): string {
  // Split path into components
  let parts = path
    .split(`.`)
    .filter((part) => ![`metrics`, `kappa_103`].includes(part) && part)

  // remove symprec value preceding rmsd
  if ([`RMSD`, `rmsd`].includes(parts.at(-1) ?? ``)) {
    parts = parts.filter((part) => !part.includes(`symprec`))
  }

  // Default formatting for other dotted paths
  return parts
    .map((part) => {
      const pretty_label = CATEGORY_LABELS[part] ?? PROPERTY_LABELS[part]
      if (pretty_label) return pretty_label
      return title_case(format_power_ten(part))
    })
    .join(` > `)
}

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
} as const // Category mapping for property paths

export const CATEGORY_LABELS: Record<string, string> = {
  discovery: `Discovery`,
  phonons: `Phonons`,
  geo_opt: `Geometry Optimization`,
  hyperparams: `Hyperparams`,
  unique_prototypes: `Unique Prototypes`,
} as const

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

export const org_logos: Record<
  string,
  { name: string; id?: string; src?: string } // Can have id OR src
> = {
  // Map of author affiliations in model YAMLs to SVG icons (either inline symbol ID
  // or external file path under /static/logos/) and full affiliation names for tooltips.
  deepmind: {
    name: `Google DeepMind`,
    src: `/logos/deepmind.svg`, // Updated path
  },
  microsoft: {
    name: `Microsoft Research`,
    id: `icon-logo-microsoft`, // Keep inline
  },
  meta: {
    name: `Meta (FAIR)`,
    id: `icon-logo-meta`, // Keep inline
  },
  cambridge: {
    name: `University of Cambridge`,
    src: `/logos/cambridge-university.svg`, // Updated path
  },
  orbital: {
    name: `Orbital Materials`,
    src: `/logos/orbital-materials.svg`, // Updated path
  },
  'seoul national university': {
    name: `Seoul National University`,
    src: `/logos/seoul-national-university.svg`, // Updated path
  },
  snu: {
    name: `Seoul National University`,
    src: `/logos/seoul-national-university.svg`,
  },
  icams: {
    name: `ICAMS, Ruhr University Bochum`,
    src: `/logos/icams-bochum.svg`, // Updated path
  },
  bochum: {
    name: `ICAMS, Ruhr University Bochum`,
    src: `/logos/icams-bochum.svg`, // Updated path
  },
  'ai for science institute': {
    name: `AI for Science Institute, Beijing`,
    src: `/logos/beijing-ai-academy.svg`, // Updated path
  },
  beijing: {
    name: `AI for Science Institute, Beijing`,
    src: `/logos/beijing-ai-academy.svg`,
  },
  cornell: {
    name: `Cornell University`,
    src: `/logos/cornell-university.svg`, // New entry
  },
  deepmd: {
    name: `DeePMD`,
    src: `/logos/deepmd.svg`, // New entry
  },
  deepmodeling: {
    name: `DeepModeling`,
    src: `/logos/deepmd.svg`, // New entry, same logo
  },
  tsinghua: {
    name: `Tsinghua University`,
    src: `/logos/tsinghua-university.svg`, // New entry
  },
  'san diego': {
    name: `UC San Diego`,
    src: `/logos/uc-san-diego.svg`, // New entry
  },
  ucsd: {
    name: `UC San Diego`,
    src: `/logos/uc-san-diego.svg`,
  },
  'massachusetts institute of technology': {
    name: `Massachusetts Institute of Technology`,
    src: `/logos/mit.svg`,
  },
  florida: {
    name: `University of Florida`,
    src: `/logos/university-of-florida.svg`,
  },
  'national institute of standards and technology': {
    name: `National Institute of Standards and Technology`,
    src: `/logos/nist.svg`,
  },
  argonne: {
    name: `Argonne National Laboratory`,
    src: `/logos/argonne-national-lab.svg`,
  },
  'university of texas at austin': {
    name: `University of Texas at Austin`,
    src: `/logos/university-of-texas-austin.svg`,
  },
  'northwestern university': {
    name: `Northwestern University`,
    src: `/logos/northwestern-university.svg`,
  },
  'chinese academy of sciences': {
    name: `Chinese Academy of Sciences`,
    src: `/logos/chinese-academy-of-sciences.svg`,
  },
  'incheon national university': {
    name: `Incheon National University`,
    src: `/logos/incheon-national-university.svg`,
  },
  'deep principle': {
    name: `Deep Principle`,
    src: `/logos/deep-principle.svg`,
  },
  'university of minnesota': {
    name: `University of Minnesota`,
    src: `/logos/university-of-minnesota.svg`,
  },
  'uc berkeley': {
    name: `University of California, Berkeley`,
    src: `/logos/uc-berkeley.svg`,
  },
}

// Attempts to find a matching logo data (ID or src) and name for a given affiliation string.
// Performs a case-insensitive search for keywords defined in org_logos.
// Returns object with logo name and either SVG ID or src path, or undefined if no match.
export function get_org_logo(
  affiliation: string,
): { name: string; id?: string; src?: string } | undefined {
  if (!affiliation) return undefined
  const lower_affiliation = affiliation.toLowerCase()
  // sort by length to prioritize longer (more specific) keys
  const sorted_keys = Object.keys(org_logos).sort((a, b) => b.length - a.length)

  for (const key of sorted_keys) {
    // Check if the lowercased affiliation string includes the lowercased key
    if (lower_affiliation.includes(key.toLowerCase())) {
      const logo_data = org_logos[key]
      return logo_data
    }
  }
}
