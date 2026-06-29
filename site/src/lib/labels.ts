import type { DiscoverySet, Label } from '$lib/types'
import MODELINGS_TASKS from '$pkg/modeling-tasks.yml'
import { ICON_DATA, type IconName } from 'matterviz'
import type {
  DatasetMetadataLabels,
  DiscoveryMetricsLabels,
  GeoOptSymmetryMetricsLabels,
  HyperparamLabels,
  MdMetricsLabels,
  MetadataLabels,
} from './schema/label.d.ts'

export const RMSD_BASELINE = 0.15 // Baseline for poor performance given worst performing model at time of writing is M3GNet at 0.1117

// Helper function to format scientific notation with superscript
// Used e.g. for symprec in geo_opt metrics
export const format_power_ten = (text: string): string =>
  text
    .replaceAll(
      /(?<base>\d+(?:\.\d+)?)e\+?(?<exponent>-?\d+)/gi,
      `$<base>×10<sup>$<exponent></sup>`,
    )
    .replace(`1×10`, `10`)

export const DISCOVERY_METRICS: DiscoveryMetricsLabels = {
  Accuracy: {
    key: `Accuracy`,
    label: `Acc`,
    description: `Accuracy of classifying crystals as thermodynamically stable or unstable`,
    better: `higher`,
    path: `metrics.discovery.unique_prototypes`,
  },
  F1: {
    key: `F1`,
    label: `F1`,
    path: `metrics.discovery.unique_prototypes`,
    description: `F1 Score: Harmonic mean of precision and recall for stable/unstable material classification`,
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
    label: `Prec`,
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
    key: `Model`,
    label: `Model`,
    description: `Model name`,
    sticky: true,
    sortable: true,
    better: undefined,
  },
  training_set: {
    key: `Training Set`,
    label: `Training Set`,
    description: `Size of and link to model training set`,
  },
  targets: {
    key: `Targets`,
    label: `Targets`,
    description: `Target property used to train the model`,
  },
  date_added: {
    key: `date_added`,
    label: `Date Added`,
    format: `%b %y`,
    description: `Submission date to the leaderboard`,
  },
  links: {
    key: `Links`,
    label: `Links`,
    description: `Model resources: paper, code repository and submission pull request`,
    sortable: false,
  },
  r_cut: {
    key: `r<sub>cut</sub>`,
    label: `r<sub>cut</sub>`,
    description: `Graph construction radius in Ångströms (cutoff distance for creating edges in the graph)`,
    unit: `Å`,
  },
  n_training_materials: {
    key: `n_training_materials`,
    label: `Training Materials`,
    description: `Number of training materials`,
    format: `~s`,
  },
  n_training_structures: {
    key: `n_training_structures`,
    label: `Training Structures`,
    description: `Number of training structures`,
    format: `~s`,
  },
  checkpoint_license: {
    key: `Ckpt License`,
    label: `Ckpt License`,
    description: `Model checkpoint license`,
    visible: false,
  },
  code_license: {
    key: `Code License`,
    label: `Code License`,
    description: `Model code license`,
    visible: false,
  },
  missing_preds: {
    key: `Missing Preds`,
    label: `Missing Preds`,
    description: `Number of missing predictions`,
    visible: false,
  },
  'Run Time (h)': {
    key: `Run Time`,
    label: `Run Time`,
    description: `Runtime in hours`,
    unit: `h`,
    visible: false,
    better: `lower`,
  },
  org: {
    key: `Org`,
    label: `Org`,
    sortable: false,
    description: `Model author affiliations`,
    cell_style: `text-align: center; width: 3.2em; max-width: 3.2em;`,
    visible: true,
    better: undefined,
  },
} as const

export const HYPERPARAMS: HyperparamLabels = {
  model_params: {
    key: `model_params`,
    label: `Params`,
    description: `Number of trainable model parameters`,
    format: `~s`,
  },
  graph_construction_radius: {
    key: `graph_construction_radius`,
    label: `r<sub>cut</sub>`,
    path: `hyperparams`,
    description: `Graph construction radius in Ångströms (cutoff distance for creating edges in the graph)`,
  },
  max_force: {
    key: `max_force`,
    label: `f<sub>max</sub>`,
    path: `hyperparams`,
    description: `Max remaining force allowed on any atom in the structure for geometry optimization convergence`,
    unit: `eV/Å`,
  },
  max_steps: {
    key: `max_steps`,
    label: `Steps`,
    path: `hyperparams`,
    description: `Maximum number of optimization steps allowed`,
  },
  ase_optimizer: {
    key: `Optimizer`,
    label: `Optimizer`,
    path: `hyperparams`,
    description: `ASE optimizer used for structure relaxation (e.g., FIRE, LBFGS, BFGS, GOQN)`,
  },
  cell_filter: {
    key: `Cell filter`,
    label: `Cell filter`,
    path: `hyperparams`,
    description: `ASE cell filter used during relaxation (e.g., FrechetCellFilter, ExpCellFilter)`,
  },
  batch_size: {
    key: `batch_size`,
    label: `Batch size`,
    path: `hyperparams`,
    description: `Batch size`,
  },
  epochs: {
    key: `epochs`,
    label: `Epochs`,
    path: `hyperparams`,
    description: `Number of training epochs`,
  },
  n_layers: {
    key: `n_layers`,
    label: `Layers`,
    path: `hyperparams`,
    description: `Number of (usually message passing) layers`,
  },
  learning_rate: {
    key: `LR`,
    label: `LR`,
    path: `hyperparams`,
    description: `Learning rate`,
  },
  max_neighbors: {
    key: `Max neighbors`,
    label: `Max neighbors`,
    path: `hyperparams`,
    description: `Maximum number of neighbors during graph construction`,
  },
  n_estimators: {
    key: `Estimators`,
    label: `Estimators`,
    path: `hyperparams`,
    description: `Number of estimators`,
  },
} as const

export const DATASET_METADATA_COLS: DatasetMetadataLabels = {
  name: { key: `Name`, label: `Name`, description: `Name of the dataset`, sticky: true },
  structures: {
    key: `Structures`,
    label: `Structures`,
    description: `Number of structures in the dataset. Any system with atomic positions and energy/force/stress labels is counted as a structure incl. successive ionic steps in MD/geometry optimization trajectories.`,
    better: `higher`,
    scale_type: `log`,
    format: `.3s`,
  },
  materials: {
    key: `Materials`,
    label: `Materials`,
    description: `Number of unique materials/prototypes in the dataset.`,
    better: `higher`,
    scale_type: `log`,
    format: `.3s`,
  },
  created: {
    key: `Created`,
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
    description: `Whether the dataset is static (fixed version) or dynamic (continuously updated).`,
    style: `text-align: center;`,
  },
  license: {
    key: `License`,
    label: `License`,
    description: `License under which the dataset is published`,
  },
  method: {
    key: `Method`,
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
    key: `Links`,
    label: `Links`,
    description: `Relevant links for the dataset`,
    sortable: false,
  },
} as const

// Object.fromEntries loses key specificity, returning Record<string, V>.
// Cast is unavoidable since the keys are dynamically constructed template literals.
export const GEO_OPT_SYMMETRY_METRICS = Object.fromEntries(
  [`1e-2`, `1e-5`]
    .flatMap((symprec) => [
      [`symmetry_match`, `=`, `higher`, `identical symmetry as`, symprec] as const,
      [`symmetry_decrease`, `↓`, `lower`, `lower symmetry than`, symprec] as const,
      [`symmetry_increase`, `↑`, null, `higher symmetry than`, symprec] as const,
    ])
    .map(([metric_key, symbol, better, desc, symprec]) => [
      `${metric_key}_${symprec}`,
      {
        key: `${metric_key}_${symprec}`,
        property: metric_key,
        symprec,
        path: `metrics.geo_opt.symprec=${symprec}`,
        label: `Σ<sub>${symbol}</sub> ${format_power_ten(symprec)}`,
        description: `Fraction of structures where ML ground state has ${desc} DFT ground state at ${format_power_ten(symprec)} symprec`,
        better,
        format: `~%`,
        // visible by default so the landing-page "Geo Opt" column preset can surface
        // them; the discovery/phonons/MD tables hide them via their own col_filter
      },
    ]),
) as unknown as GeoOptSymmetryMetricsLabels

export const MD_METRICS: MdMetricsLabels = {
  md_energy_rmse: {
    key: `energy_rmse`,
    label: `ΔE<sub>RMSE</sub>`,
    description: `Root mean squared error of model vs ab-initio energy fluctuations (each trajectory's mean energy is subtracted from both) on reference MD frames. Mean-subtraction makes this invariant to the absolute energy reference, which differs across CFPMD-26 systems (all-electron vs PAW)`,
    unit: `meV/atom`,
    path: `metrics.md`,
    better: `lower`,
    format: `.1f`,
    style: `border-left: 1px solid black;`,
  },
  md_force_rmse: {
    key: `force_rmse`,
    label: `F<sub>RMSE</sub>`,
    description: `Root mean squared error of model-predicted forces on reference ab-initio MD frames`,
    unit: `meV/Å`,
    path: `metrics.md`,
    better: `lower`,
    format: `.1f`,
  },
  md_rdf_error: {
    key: `rdf_error`,
    label: `ΔRDF`,
    description: `Mean radial distribution function error between MLIP and ab-initio MD trajectories. 0% = perfect match, 100% = as different from the reference as an ideal gas`,
    unit: `%`,
    path: `metrics.md`,
    range: [0, 100],
    better: `lower`,
    format: `.1f`,
  },
  md_adf_error: {
    key: `adf_error`,
    label: `ΔADF`,
    description: `Bond-angle distribution error between MLIP and ab-initio MD trajectories. Wasserstein-1 distance over first-coordination-shell bond angles (neighbor pairs within species-aware covalent-radius cutoffs), normalized by the reference's distance to a structureless background. 0% = perfect match, 100% = as different from the reference as a featureless angular distribution`,
    unit: `%`,
    path: `metrics.md`,
    range: [0, 100],
    better: `lower`,
    format: `.1f`,
  },
  md_vdos_error: {
    key: `vdos_error`,
    label: `ΔvDOS`,
    description: `Vibrational density of states error between MLIP and ab-initio MD trajectories (from the velocity autocorrelation spectrum): Wasserstein-1 distance between the spectra normalized by the reference's spectral spread, so it grows with systematic frequency softening/hardening rather than saturating. 0% = perfect match, 100% = displaced by at least the reference's own spectral width`,
    unit: `%`,
    path: `metrics.md`,
    range: [0, 100],
    better: `lower`,
    format: `.1f`,
  },
  md_pressure_mae: {
    key: `pressure_mae`,
    label: `P<sub>MAE</sub>`,
    description: `Absolute difference between the mean pressures of the MLIP and ab-initio MD trajectories (mean-stress bias); independent of frame pairing`,
    unit: `GPa`,
    path: `metrics.md`,
    better: `lower`,
    format: `.2f`,
  },
  md_pressure_wasserstein: {
    key: `pressure_wasserstein`,
    label: `P<sub>W1</sub>`,
    description: `Wasserstein-1 distance between pressure distributions of MLIP and ab-initio MD trajectories (insensitive to frame pairing)`,
    unit: `GPa`,
    path: `metrics.md`,
    better: `lower`,
    format: `.2f`,
  },
  md_pressure_error: {
    key: `pressure_error`,
    label: `ΔP`,
    description: `Pressure-distribution error: the non-overlap of the area-normalized MLIP and ab-initio pressure histograms over shared bin edges. 0% = identical distributions, 100% = disjoint`,
    unit: `%`,
    path: `metrics.md`,
    range: [0, 100],
    better: `lower`,
    format: `.1f`,
  },
  md_combined_score: {
    key: `combined_score`,
    label: `CMDS`,
    description: `Combined MD score in [0,1] (higher is better): 1 − mean(ΔRDF, ΔADF, ΔvDOS, ΔP) / 100, where each error is a percentage and lower is better. Higher = closer to ab-initio dynamics; intended to feed into CPS as a normalized component.`,
    path: `metrics.md`,
    range: [0, 1],
    better: `higher`,
    format: `.3f`,
  },
} as const

type DiatomicsMetricKey =
  | `tortuosity`
  | `energy_diff_flips`
  | `energy_grad_norm_max`
  | `energy_jump`
  | `conservation`
  | `force_flips`
  | `force_total_variation`
  | `force_jump`

export const DIATOMICS_METRICS: Record<DiatomicsMetricKey, Label> = {
  tortuosity: {
    key: `tortuosity`,
    label: `τ`,
    description: `Mean tortuosity of homonuclear diatomic energy curves, measuring extra energy variation beyond a single-well or monotonic curve`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
    style: `border-left: 1px solid black;`,
  },
  energy_diff_flips: {
    key: `energy_diff_flips`,
    label: `E flips`,
    description: `Mean number of sign flips in adjacent diatomic energy differences`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
  },
  energy_grad_norm_max: {
    key: `energy_grad_norm_max`,
    label: `max |∇E|`,
    description: `Mean maximum absolute energy gradient along homonuclear diatomic curves`,
    unit: `eV/Å`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
  },
  energy_jump: {
    key: `energy_jump`,
    label: `E jump`,
    description: `Mean energy jump at sign-flip points in homonuclear diatomic curves`,
    unit: `eV`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
  },
  conservation: {
    key: `conservation`,
    label: `Conserv.`,
    description: `Mean deviation between predicted forces and the negative gradient of predicted diatomic energies`,
    unit: `eV/Å`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
  },
  force_flips: {
    key: `force_flips`,
    label: `F flips`,
    description: `Mean number of force-direction flips along homonuclear diatomic curves`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
  },
  force_total_variation: {
    key: `force_total_variation`,
    label: `F TV`,
    description: `Mean total variation of forces along homonuclear diatomic curves`,
    unit: `eV/Å`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
  },
  force_jump: {
    key: `force_jump`,
    label: `F jump`,
    description: `Mean force jump at force-direction flip points in homonuclear diatomic curves`,
    unit: `eV/Å`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
  },
}

export type AllMetrics = DiscoveryMetricsLabels &
  GeoOptSymmetryMetricsLabels &
  MdMetricsLabels &
  Record<DiatomicsMetricKey, Label> & {
    CPS: Label
    κ_SRME: Label
    κ_SRE: Label
    RMSD: Label
  }

export const ALL_METRICS: AllMetrics = {
  // Dynamic metrics
  CPS: {
    key: `CPS`,
    label: `CPS`,
    description: `Combined Performance Score averages discovery (F1), structure optimization (RMSD), and phonon performance (κ<sub>SRME</sub>) according to user-defined weights. Warning: This is not a stable metric. Further prediction tasks will be added to it in the future with the goal of making it a more holistic measure of overall model utility over time. When referring to it in papers, best include the benchmark version to avoid confusion (e.g. CPS-1 for the first version of CPS introduced in Matbench Discovery v1)`,
    range: [0, 1],
    better: `higher`,
    format: `.3f`,
  },
  ...DISCOVERY_METRICS,
  // Phonon metrics
  κ_SRME: {
    key: `κ_SRME`,
    label: `κ<sub>SRME</sub>`,
    description: `Symmetric relative mean error in predicted phonon mode contributions to thermal conductivity κ`,
    path: `metrics.phonons.kappa_103`,
    better: `lower`,
  },
  κ_SRE: {
    key: `κ_SRE`,
    label: `κ<sub>SRE</sub>`,
    description: `Symmetric relative error of total lattice thermal conductivity κ, averaged over the 103 phononDB-PBE materials (range [0, 2])`,
    path: `metrics.phonons.kappa_103`,
    better: `lower`,
  },
  // Geometry optimization metrics
  RMSD: {
    key: `rmsd`,
    path: `metrics.geo_opt.symprec=1e-2`,
    label: `RMSD`,
    // Unit intentionally hidden for concise table column headers.
    range: [0, RMSD_BASELINE],
    better: `lower`,
    description: `Normalized, unitless StructureMatcher RMSD between ML- and DFT-relaxed structures after matching; unmatched structures are assigned 1.0`,
    style: `border-left: 1px solid black;`,
  },
  ...GEO_OPT_SYMMETRY_METRICS,
  ...MD_METRICS,
  ...DIATOMICS_METRICS,
} as const

export const DISCOVERY_SET_LABELS: Record<
  DiscoverySet,
  { label: string; description: string; link?: string }
> = {
  full_test_set: {
    label: `Full Test Set`,
    description: `<strong>257k total structures</strong><br/>
      Metrics computed on all WBM structures including duplicate structure prototypes`,
  },
  unique_prototypes: {
    label: `Unique Prototypes`,
    description: `<strong>~215k unique prototypes</strong><br/>
      Deduplicated by matching Aflow-style prototypes.<br/>
      Use this to avoid counting similar structures that should relax to same ground state multiple times.`,
    link: `https://github.com/janosh/matbench-discovery/blob/37baf7986f848/data/wbm/compile_wbm_test_set.py#L640-L654`,
  },
  most_stable_10k: {
    label: `10k Most Stable`,
    description: `<strong>Top 10k predictions by model</strong><br/>
      Each model's structures by lowest predicted energy above hull.<br/>
      Use this to evaluate discovery performance in an actual discovery campaign at fixed compute budget e.g. for DFT validation.`,
  },
} as const

// SelectToggle options for switching between WBM test subsets
export const discovery_set_toggle_options = Object.entries(DISCOVERY_SET_LABELS).map(
  ([value, { label, description: tooltip, link }]) => ({ value, label, tooltip, link }),
)

const PROPERTY_LABELS = Object.fromEntries(
  Object.values({ ...ALL_METRICS, ...METADATA_COLS, ...HYPERPARAMS }).map((prop) => [
    prop.key ?? prop.label,
    prop.label,
  ]),
)

// Formats a property path for display in UI components
export function format_property_path(path: string): string {
  // Split path into components
  let parts = path
    .split(`.`)
    .filter((part) => ![`metrics`, `kappa_103`].includes(part) && part)

  // Remove symprec value preceding rmsd
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

const CATEGORY_LABELS = Object.fromEntries(
  Object.entries({ ...MODELINGS_TASKS, ...DISCOVERY_SET_LABELS, ...HYPERPARAMS }).map(
    ([key, task]) => [key, task.label],
  ),
)
// Add explicit mapping for hyperparams to show as "Hyperparams"
CATEGORY_LABELS.hyperparams = `Hyperparams`

const to_title = (str: string) => str.charAt(0).toUpperCase() + str.slice(1)
export const title_case = (str: string) =>
  str.replaceAll(`_`, ` `).split(` `).map(to_title).join(` `)

// Return type for get_org_logo function
export interface OrgLogo {
  name: string
  id?: string
  src?: string
  validated_icon?: IconName
}

// Map of author affiliations in model YAMLs to SVG icons (either inline symbol ID
// Or external file path under /static/logos/) and full affiliation names for tooltips.
const org_logos = {
  'AI for Science Institute, Beijing': `/logos/beijing-ai-for-science-institute.svg`,
  'Argonne National Laboratory': `/logos/argonne-national-lab.svg`,
  'Beijing Institute of Applied Physics and Computational Mathematics (IAPCM)': `/logos/beijing-iapcm.svg`,
  'Chinese Academy of Sciences': `/logos/chinese-academy-of-sciences.svg`,
  'Cornell University': `/logos/cornell-university.svg`,
  'DAMO Academy, Alibaba Inc': `/logos/damo-alibaba-logo.svg`,
  'Deep Principle': `/logos/deep-principle.svg`,
  DeePMD: `/logos/deepmd.svg`,
  'FAIR at Meta': `icon:LogoMeta`,
  'Google DeepMind': `/logos/deepmind.svg`,
  'ICAMS, Ruhr University Bochum': `/logos/icams-bochum.svg`,
  'Incheon National University': `/logos/incheon-national-university.svg`,
  'Institute of Computing Technology, Chinese Academy of Science, Beijing': `/logos/ict-cas-beijing.svg`,
  'Massachusetts Institute of Technology': `/logos/mit.svg`,
  'Microsoft Research': `icon:LogoMicrosoft`,
  'MIR Group, Harvard University': `/logos/mir-group-harvard.svg`,
  'Mirror Physics': `/logos/mirror-physics.svg`,
  'National Institute of Standards and Technology': `/logos/nist.svg`,
  'Ningbo Institute of Artificial Intelligence Industry': `/logos/ningbo-institute-of-artificial-intelligence-industry.svg`,
  'Northwestern University': `/logos/northwestern-university.svg`,
  'Orbital Materials': `/logos/orbital-materials.svg`,
  'Peking University': `/logos/peking-university.svg`,
  'Seoul National University': `/logos/seoul-national-university.svg`,
  'Samsung Electronics': `/logos/samsung-electronics.svg`,
  'Texas A&M University': `/logos/texas-a&m.svg`,
  'Tsinghua University': `/logos/tsinghua-university.svg`,
  'UC San Diego': `/logos/uc-san-diego.svg`,
  'UC Berkeley': `/logos/uc-berkeley.svg`,
  'University of California, Los Angeles': `/logos/ucla.svg`,
  'University of Cambridge': `/logos/cambridge-university.svg`,
  'University of Florida': `/logos/university-of-florida.svg`,
  'University of Minnesota': `/logos/university-of-minnesota.svg`,
  'University of Texas at Austin': `/logos/university-of-texas-austin.svg`,
  'Beijing Information Science and Technology University': `/logos/beijing-information-science-and-technology-university.svg`,
  'Zhejiang Lab': `/logos/zhejiang-lab.svg`,
  EPFL: `/logos/epfl.svg`,
  'ShanghaiTech University': `/logos/shanghaitech-university.svg`,
  'Nanjing University': `/logos/nanjing-university.svg`,
} as const

// Attempts to find a matching logo data (ID or src) and name for a given affiliation string.
// Performs a case-insensitive search for keywords defined in org_logos.
// Returns object with logo name and either SVG ID or src path, or undefined if no match.
// For icon references (icon:*), validates against matterviz IconName and includes validated_icon.
export function get_org_logo(affiliation: string): OrgLogo | undefined {
  if (!affiliation) return undefined

  for (const [key_val, logo_val] of Object.entries(org_logos)) {
    // Check if lowercased affiliation string includes lowercased org key
    if (affiliation.toLowerCase().includes(key_val.toLowerCase())) {
      if (logo_val.startsWith(`/logos/`)) {
        return { name: key_val, src: logo_val }
      } else if (logo_val.startsWith(`icon:`)) {
        const icon_name = logo_val.replace(`icon:`, ``)
        const validated_icon =
          icon_name in ICON_DATA ? (icon_name as IconName) : undefined
        if (!validated_icon && !import.meta.env.PROD) {
          console.warn(
            `Invalid icon name "${icon_name}" for org "${key_val}". ` +
              `Valid names: ${Object.keys(ICON_DATA).slice(0, 10).join(`, `)}...`,
          )
        }
        return { name: key_val, id: logo_val, validated_icon }
      }
      return { name: key_val, id: logo_val }
    }
  }
  return undefined
}
