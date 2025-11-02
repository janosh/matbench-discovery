import type { DiscoverySet, Label } from '$lib/types'
import MODELINGS_TASKS from '$pkg/modeling-tasks.yml'
import type {
  DatasetMetadataLabels,
  DiscoveryMetricsLabels,
  GeoOptSymmetryMetricsLabels,
  HyperparamLabels,
  MetadataLabels,
} from './label-schema.d.ts'

export const RMSD_BASELINE = 0.15 // baseline for poor performance given worst performing model at time of writing is M3GNet at 0.1117

// Helper function to format scientific notation with superscript
// used e.g. for symprec in geo_opt metrics
export const format_power_ten = (text: string): string =>
  text
    .replace(/(\d+(?:\.\d+)?)e\+?(-?\d+)/gi, (_, base, exponent) =>
      `${base}×10<sup>${exponent}</sup>`)
    .replace(`1×10`, `10`)

export const DISCOVERY_METRICS: DiscoveryMetricsLabels = {
  Accuracy: {
    key: `Accuracy`,
    short: `Acc`,
    label: `Accuracy`,
    description:
      `Accuracy of classifying crystals as thermodynamically stable or unstable`,
    better: `higher`,
    path: `metrics.discovery.unique_prototypes`,
  },
  F1: {
    key: `F1`,
    short: `F1`,
    label: `F1 Score`,
    path: `metrics.discovery.unique_prototypes`,
    description:
      `Harmonic mean of precision and recall for stable/unstable material classification`,
    range: [0, 1],
    better: `higher`,
  },
  DAF: {
    key: `DAF`,
    label: `DAF`,
    description:
      `Discovery Acceleration Factor measuring how much better ML models classify thermodynamic stability compared to random guessing`,
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
    description: `Mean Absolute Error of predicted vs. DFT convex hull distance`,
    unit: `eV / atom`,
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
    unit: `eV / atom`,
    path: `metrics.discovery.unique_prototypes`,
    better: `lower`,
  },
} as const

export const METADATA_COLS: MetadataLabels = {
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
    description:
      `Graph construction radius in Ångströms (cutoff distance for creating edges in the graph)`,
    unit: `Å`,
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
  'Run Time (h)': {
    key: `run_time_h`,
    label: `Run Time`,
    description: `Runtime in hours`,
    unit: `h`,
    visible: false,
    better: `lower`,
  },
  org: {
    key: `org`,
    label: `Org`,
    sortable: false,
    description: `Most common author affiliations`,
    cell_style: `text-align: center; max-width: 2em; transform: scale(1.3);`,
    visible: true,
    better: null,
  },
} as const

export const HYPERPARAMS: HyperparamLabels = {
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
    description:
      `Graph construction radius in Ångströms (cutoff distance for creating edges in the graph)`,
  },
  max_force: {
    label: `Max force`,
    key: `max_force`,
    short: `f<sub>max</sub>`,
    path: `hyperparams`,
    description:
      `Max remaining force allowed on any atom in the structure for geometry optimization convergence`,
    unit: `eV/Å`,
  },
  max_steps: {
    label: `Max relaxation steps`,
    key: `max_steps`,
    short: `Steps`,
    path: `hyperparams`,
    description: `Maximum number of optimization steps allowed`,
  },
  ase_optimizer: {
    label: `ASE optimizer`,
    key: `ase_optimizer`,
    short: `Optimizer`,
    path: `hyperparams`,
    description:
      `ASE optimizer used for structure relaxation (e.g., FIRE, LBFGS, BFGS, GOQN)`,
  },
  cell_filter: {
    label: `Cell filter`,
    key: `cell_filter`,
    path: `hyperparams`,
    description:
      `ASE cell filter used during relaxation (e.g., FrechetCellFilter, ExpCellFilter)`,
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

export const DATASET_METADATA_COLS: DatasetMetadataLabels = {
  name: { key: `name`, label: `Name`, description: `Name of the dataset`, sticky: true },
  structures: {
    key: `n_structures`,
    label: `Number of Structures`,
    short: `Structures`,
    description:
      `Number of structures in the dataset. Any system with atomic positions and energy/force/stress labels is counted as a structure incl. successive ionic steps in MD/geometry optimization trajectories.`,
    better: `higher`,
    scale_type: `log`,
    format: `.3s`,
  },
  materials: {
    key: `n_materials`,
    label: `Number of Materials`,
    short: `Materials`,
    description: `Number of unique materials/prototypes in the dataset.`,
    better: `higher`,
    scale_type: `log`,
    format: `.3s`,
  },
  created: {
    key: `created`,
    label: `Created`,
    description: `Date the dataset was created/started`,
  },
  open: {
    key: `Open`,
    label: `Open`,
    description: `Whether the dataset is openly available`,
    style: `text-align: center;`,
  },
  static: {
    key: `Static`,
    label: `Static`,
    description:
      `Whether the dataset is static (fixed version) or dynamic (continuously updated).`,
    style: `text-align: center;`,
  },
  license: {
    key: `license`,
    label: `License`,
    description: `License under which the dataset is published`,
  },
  method: {
    key: `method`,
    label: `Method`,
    description: `Method(s) used to generate the data`,
    style: `max-width: 5em;`,
  },
  api: {
    key: `API`,
    label: `API`,
    description: `API docs (OPTIMADE or native)`,
    sortable: false,
  },
  links: {
    key: `links`,
    label: `Links`,
    description: `Relevant links for the dataset`,
    sortable: false,
  },
} as const

export const GEO_OPT_SYMMETRY_METRICS = Object.fromEntries(
  [`1e-2`, `1e-5`]
    .flatMap(
      (symprec) => [
        [`symmetry_match`, `=`, `higher`, `identical symmetry as`, symprec] as const,
        [`symmetry_decrease`, `↓`, `lower`, `lower symmetry than`, symprec] as const,
        [`symmetry_increase`, `↑`, null, `higher symmetry than`, symprec] as const,
      ],
    )
    .map(([key, symbol, better, desc, symprec]) => [
      `${key}_${symprec}`,
      {
        key,
        symprec,
        path: `metrics.geo_opt.symprec=${symprec}`,
        short: `Σ<sub>${symbol}</sub> ${format_power_ten(symprec)}`,
        label: `Σ<sub>${symbol}</sub> (symprec=${format_power_ten(symprec)})`,
        description:
          `Fraction of structures where ML ground state has ${desc} DFT ground state at ${
            format_power_ten(symprec)
          } symprec`,
        better,
        format: `~%`,
        visible: false,
      },
    ]),
) as unknown as GeoOptSymmetryMetricsLabels

export type AllMetrics =
  & DiscoveryMetricsLabels
  & GeoOptSymmetryMetricsLabels
  & { CPS: Label; κ_SRME: Label; RMSD: Label }

export const ALL_METRICS: AllMetrics = {
  // Dynamic metrics
  CPS: {
    key: `CPS`,
    short: `CPS`,
    label: `Combined Performance Score`,
    description:
      `Combined Performance Score averages discovery (F1), structure optimization (RMSD), and phonon performance (κ<sub>SRME</sub>) according to user-defined weights. Warning: This is not a stable metric. Further prediction tasks will be added to it in the future with the goal of making it a more holistic measure of overall model utility over time. When referring to it in papers, best include the benchmark version to avoid confusion (e.g. CPS-1 for the first version of CPS introduced in Matbench Discovery v1)`,
    range: [0, 1],
    better: `higher`,
    format: `.3f`,
  },
  ...DISCOVERY_METRICS,
  // Phonon metrics
  κ_SRME: {
    key: `κ_SRME`,
    label: `κ<sub>SRME</sub>`,
    description:
      `Symmetric relative mean error in predicted phonon mode contributions to thermal conductivity κ`,
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
    description:
      `Root mean squared displacement between predicted and reference structures after relaxation`,
    style: `border-left: 1px solid black;`,
  },
  ...GEO_OPT_SYMMETRY_METRICS,
} as const

export const DISCOVERY_SET_LABELS: Record<
  DiscoverySet,
  { label: string; description: string; link?: string }
> = {
  full_test_set: {
    label: `Full Test Set`,
    description:
      `Metrics computed on the full 257k WBM test set including duplicate structure prototypes`,
  },
  unique_prototypes: {
    label: `Unique Prototypes`,
    description:
      `Metrics computed only on ~215k unique structure prototypes in WBM determined by matching Aflow-style prototype strings.`,
    link:
      `https://github.com/janosh/matbench-discovery/blob/37baf7986f848/data/wbm/compile_wbm_test_set.py#L640-L654`,
  },
  most_stable_10k: {
    label: `10k Most Stable`,
    description:
      `Metrics computed on the 10k structures predicted to be most stable (different for each model)`,
  },
} as const

export const PROPERTY_LABELS = Object.fromEntries(
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
// Add explicit mapping for hyperparams to show as "Hyperparams"
CATEGORY_LABELS.hyperparams = `Hyperparams`

export const to_title = (str: string) => str.charAt(0).toUpperCase() + str.slice(1)
export const title_case = (str: string) =>
  str.replaceAll(`_`, ` `).split(` `).map(to_title).join(` `)

// Map of author affiliations in model YAMLs to SVG icons (either inline symbol ID
// or external file path under /static/logos/) and full affiliation names for tooltips.
export const org_logos = {
  'AI for Science Institute, Beijing': `/logos/beijing-ai-for-science-institute.svg`,
  'Argonne National Laboratory': `/logos/argonne-national-lab.svg`,
  'Beijing Institute of Applied Physics and Computational Mathematics (IAPCM)':
    `/logos/beijing-iapcm.svg`,
  'Chinese Academy of Sciences': `/logos/chinese-academy-of-sciences.svg`,
  'Cornell University': `/logos/cornell-university.svg`,
  'DAMO Academy, Alibaba Inc': `/logos/damo-alibaba-logo.svg`,
  'Deep Principle': `/logos/deep-principle.svg`,
  DeePMD: `/logos/deepmd.svg`,
  'FAIR at Meta': `icon:Meta`,
  'Google DeepMind': `/logos/deepmind.svg`,
  'ICAMS, Ruhr University Bochum': `/logos/icams-bochum.svg`,
  'Incheon National University': `/logos/incheon-national-university.svg`,
  'Institute of Computing Technology, Chinese Academy of Science, Beijing':
    `/logos/ict-cas-beijing.svg`,
  'Massachusetts Institute of Technology': `/logos/mit.svg`,
  'Microsoft Research': `icon:Microsoft`,
  'MIR Group, Harvard University': `/logos/mir-group-harvard.svg`,
  'National Institute of Standards and Technology': `/logos/nist.svg`,
  'Ningbo Institute of Artificial Intelligence Industry':
    `/logos/ningbo-institute-of-artificial-intelligence-industry.svg`,
  'Northwestern University': `/logos/northwestern-university.svg`,
  'Orbital Materials': `/logos/orbital-materials.svg`,
  'Seoul National University': `/logos/seoul-national-university.svg`,
  'Materials AI Lab at Samsung Electronics': `/logos/samsung-electronics.svg`,
  'Texas A&M University': `/logos/texas-a&m.svg`,
  'Tsinghua University': `/logos/tsinghua-university.svg`,
  'UC San Diego': `/logos/uc-san-diego.svg`,
  'UC Berkeley': `/logos/uc-berkeley.svg`,
  'University of Cambridge': `/logos/cambridge-university.svg`,
  'University of Florida': `/logos/university-of-florida.svg`,
  'University of Minnesota': `/logos/university-of-minnesota.svg`,
  'University of Texas at Austin': `/logos/university-of-texas-austin.svg`,
  'Beijing Information Science and Technology University':
    `/logos/beijing-information-science-and-technology-university.svg`,
  'Zhejiang Lab': `/logos/zhejiang-lab.svg`,
} as const

// Attempts to find a matching logo data (ID or src) and name for a given affiliation string.
// Performs a case-insensitive search for keywords defined in org_logos.
// Returns object with logo name and either SVG ID or src path, or undefined if no match.
export function get_org_logo(
  affiliation: string,
): { name: string; id?: string; src?: string } | undefined {
  if (!affiliation) return undefined

  for (const [key_val, logo_val] of Object.entries(org_logos)) {
    // Check if lowercased affiliation string includes lowercased org key
    if (affiliation.includes(key_val)) {
      if (logo_val.startsWith(`/logos/`)) return { name: key_val, src: logo_val }
      else return { name: key_val, id: logo_val }
    }
  }
  return undefined
}
