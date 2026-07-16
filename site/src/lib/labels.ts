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
  training_sets: {
    key: `Training Set`,
    label: `Training Set`,
    description: `Size of and link to model training set`,
  },
  targets: {
    key: `Targets`,
    label: `Targets`,
    description: `Target property used to train the model`,
  },
  benchmark_added: {
    key: `benchmark_added`,
    label: `Date Added`,
    path: `dates`,
    format: `%b %y`,
    description: `Date the model was included on the benchmark leaderboard`,
  },
  links: {
    key: `Links`,
    label: `Links`,
    description: `Model resources: paper, code repository and submission pull request`,
    sortable: false,
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
    cell_style: `text-align: center; width: 2.5em; max-width: 2.5em;`,
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
    path: `hyperparams.architecture`,
    description: `Graph construction radius in Ångströms (cutoff distance for creating edges in the graph)`,
  },
  max_force: {
    key: `max_force`,
    label: `f<sub>max</sub>`,
    path: `hyperparams.evaluation`,
    description: `Max remaining force allowed on any atom in the structure for geometry optimization convergence`,
    unit: `eV/Å`,
  },
  max_steps: {
    key: `max_steps`,
    label: `Steps`,
    path: `hyperparams.evaluation`,
    description: `Maximum number of optimization steps allowed`,
  },
  ase_optimizer: {
    key: `Optimizer`,
    label: `Optimizer`,
    path: `hyperparams.evaluation`,
    description: `ASE optimizer used for structure relaxation (e.g., FIRE, LBFGS, BFGS, GOQN)`,
  },
  cell_filter: {
    key: `Cell filter`,
    label: `Cell filter`,
    path: `hyperparams.evaluation`,
    description: `ASE cell filter used during relaxation (e.g., FrechetCellFilter, ExpCellFilter)`,
  },
  batch_size: {
    key: `batch_size`,
    label: `Batch size`,
    path: `hyperparams.training`,
    description: `Batch size`,
  },
  epochs: {
    key: `epochs`,
    label: `Epochs`,
    path: `hyperparams.training`,
    description: `Number of training epochs`,
  },
  n_layers: {
    key: `n_layers`,
    label: `Layers`,
    path: `hyperparams.architecture`,
    description: `Number of (usually message passing) layers`,
  },
  learning_rate: {
    key: `LR`,
    label: `LR`,
    path: `hyperparams.training`,
    description: `Learning rate`,
  },
  max_neighbors: {
    key: `Max neighbors`,
    label: `Max neighbors`,
    path: `hyperparams.architecture`,
    description: `Maximum number of neighbors during graph construction`,
  },
  n_estimators: {
    key: `n_estimators`,
    label: `Estimators`,
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
    description: `Private-label diagnostic: root mean squared error of model vs ab-initio energy fluctuations on reference MD frames (each trajectory's mean energy is subtracted from both). Excluded from CMDS because public DynaMat reference files do not include energies`,
    unit: `meV/atom`,
    path: `metrics.md`,
    better: `lower`,
    format: `.1f`,
    style: `border-left: 1px solid black;`,
  },
  md_force_rmse: {
    key: `force_rmse`,
    label: `F<sub>RMSE</sub>`,
    description: `Private-label diagnostic: root mean squared error of model-predicted forces on reference ab-initio MD frames. Excluded from CMDS because public DynaMat reference files do not include forces`,
    unit: `meV/Å`,
    path: `metrics.md`,
    better: `lower`,
    format: `.1f`,
  },
  md_rdf_error: {
    key: `rdf_error`,
    label: `ΔRDF`,
    description: `Mean radial distribution function error between MLIP and ab-initio MD trajectories. 0% = perfect match, 100% = as different from the reference as an ideal gas. Hidden from the leaderboard and excluded from CMDS: it correlates 0.9+ with ΔvDOS and ΔADF across models, so it adds columns but little signal`,
    unit: `%`,
    path: `metrics.md`,
    range: [0, 100],
    better: `lower`,
    format: `.1f`,
    visible: false,
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
    description: `Combined MD score in [0,1] (higher is better): weighted mean of the ΔvDOS (30%), ΔADF (20%), ΔP (30%) subscores (1 − error/100) and Speed (summed rollout wall time, log-scaled, 20%), reweightable on the MD task page. Computed on the fly like CPS/CDS, never stored with submissions. ΔRDF is excluded as redundant (0.9+ correlation with ΔvDOS/ΔADF); models without recorded timings get no CMDS unless the Speed weight is zeroed. Higher = closer to ab-initio dynamics.`,
    path: `metrics.md`,
    range: [0, 1],
    better: `higher`,
    format: `.3f`,
  },
  md_run_time_sec: {
    key: `md_run_time_sec`,
    property: `run_time_sec`,
    label: `Speed`,
    description: `MD wall time in seconds to roll out all 17 DynaMat v1.0 NVT trajectories (20 ps each), summed over systems, excluding metric evaluation. All timings to date were measured on a single NVIDIA H200 per system (recorded in the model YAML's hardware field); blank for submissions without recorded timings`,
    unit: `s`,
    path: `metrics.md`,
    better: `lower`,
    format: `.3~s`,
  },
  md_time_multiplier: {
    key: `md_time_multiplier`,
    label: `Slowdown`,
    description: `MD wall time as a multiple of the fastest finite MD wall time among models matching the active task, energy-only, training-data and openness filters (1× = fastest model in view)`,
    unit: `×`,
    better: `lower`,
    format: `.2~f`,
  },
  // memory columns are hidden by default until enough models have re-run with memory
  // tracking (published all-or-nothing over the 17 systems); toggle via column controls
  md_max_gpu_mem_gb: {
    key: `md_max_gpu_mem_gb`,
    property: `max_gpu_mem_gb`,
    label: `VRAM`,
    description: `Peak GPU memory (torch CUDA allocator high-water mark) over all 17 DynaMat v1.0 NVT rollouts, set by the largest system (~500 atoms). Answers "what GPU does this model need for MD?"`,
    unit: `GB`,
    path: `metrics.md`,
    better: `lower`,
    format: `.3~s`,
    visible: false,
  },
  md_max_rss_gb: {
    key: `md_max_rss_gb`,
    property: `max_rss_gb`,
    label: `RAM`,
    description: `Peak host memory (resident set size high-water mark) over all 17 DynaMat v1.0 NVT rollouts`,
    unit: `GB`,
    path: `metrics.md`,
    better: `lower`,
    format: `.3~s`,
    visible: false,
  },
} as const

type DiatomicsMetricKey =
  | `tortuosity`
  | `energy_diff_flips`
  | `energy_jump`
  | `pbe_wall_dist_mae`
  | `pbe_energy_mae`
  | `pbe_bond_length_error`
  | `pbe_well_depth_error`
  | `pbe_force_mae`
  | `pbe_vib_freq_error`
  | `force_flips`
  | `force_total_variation`
  | `force_jump`
  | `diatomics_combined_score`
  | `diatomics_run_time_sec`
  | `diatomics_time_multiplier`
  | `diatomics_max_gpu_mem_gb`
  | `diatomics_max_rss_gb`

const scored_diatomic_range = `scored range from 0.9× covalent radius to min(3.1× Alvarez vdW radius, max sampled distance)`
const scored_diatomic_wall_range = `repulsive-wall range from 0.8× covalent radius to min(3.1× Alvarez vdW radius, max sampled distance)`

export const DIATOMICS_METRICS: Record<DiatomicsMetricKey, Label> = {
  tortuosity: {
    key: `tortuosity`,
    label: `τ`,
    description: `Mean tortuosity over ${scored_diatomic_range}, measuring extra energy variation beyond a single-well or monotonic curve`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
    style: `border-left: 1px solid black;`,
  },
  energy_diff_flips: {
    key: `energy_diff_flips`,
    label: `E flips`,
    description: `Mean number of sign flips in adjacent diatomic energy differences over ${scored_diatomic_range}`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
  },
  energy_jump: {
    key: `energy_jump`,
    label: `E jump`,
    description: `Mean energy jump at sign-flip points over ${scored_diatomic_range}`,
    unit: `eV`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
  },
  pbe_wall_dist_mae: {
    key: `pbe_wall_dist_mae`,
    label: `PBE Δr wall`,
    description: `Mean repulsive-wall distance error relative to PBE over the ${scored_diatomic_wall_range}, at 1, 5, 10, 20, 50 and 100 eV above the well minimum where reached by the reference; predictions that miss a supported threshold receive the full reference-radius error`,
    unit: `Å`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
  },
  pbe_energy_mae: {
    key: `pbe_energy_mae`,
    label: `PBE E MAE`,
    description: `Mean absolute energy error relative to PBE over ${scored_diatomic_range}, after aligning curves at the largest shared separation`,
    unit: `eV`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
  },
  pbe_bond_length_error: {
    key: `pbe_bond_length_error`,
    label: `PBE Δr<sub>e</sub>`,
    description: `Equilibrium bond-length error relative to PBE from a local quadratic fit near the minimum over ${scored_diatomic_range}`,
    unit: `Å`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
  },
  pbe_well_depth_error: {
    key: `pbe_well_depth_error`,
    label: `PBE ΔD<sub>e</sub>`,
    description: `Well-depth error relative to PBE over ${scored_diatomic_range}, D_e = E(r_max) - E_min`,
    unit: `eV`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
  },
  pbe_force_mae: {
    key: `pbe_force_mae`,
    label: `PBE F MAE`,
    description: `Mean absolute force error relative to PBE forces over ${scored_diatomic_range}`,
    unit: `eV/Å`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
  },
  pbe_vib_freq_error: {
    key: `pbe_vib_freq_error`,
    label: `PBE Δω`,
    description: `Harmonic vibrational-frequency error relative to PBE from a local quadratic fit near the minimum over ${scored_diatomic_range}`,
    unit: `cm⁻¹`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
  },
  force_flips: {
    key: `force_flips`,
    label: `F flips`,
    description: `Mean number of force-direction flips over ${scored_diatomic_range}`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
  },
  force_total_variation: {
    key: `force_total_variation`,
    label: `F TV`,
    description: `Mean total variation of forces over ${scored_diatomic_range}`,
    unit: `eV/Å`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
  },
  force_jump: {
    key: `force_jump`,
    label: `F jump`,
    description: `Mean force jump at force-direction flip points over ${scored_diatomic_range}`,
    unit: `eV/Å`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~g`,
  },
  diatomics_combined_score: {
    key: `diatomics_combined_score`,
    property: `combined_score`,
    label: `CDS`,
    description: `Combined Diatomics Score in [0,1] (higher is better): weighted mean of 4 pillar subscores - Accuracy (PBE energy/force MAE, 44%), Geometry (repulsive wall, bond length, well depth vs PBE, 22%), Physicality (reference-free smoothness: energy jumps, force flips, tortuosity, 22%) and Speed (sweep wall time, log-scaled, 11%; hardware varies by submission) - multiplied by the fraction of 87 benchmark elements completed and reweightable on the diatomics page. Components use fixed normalization baselines so scores remain stable as models are added. Computed on the fly like CPS/CMDS, never stored with submissions`,
    path: `metrics.diatomics`,
    range: [0, 1],
    better: `higher`,
    format: `.3f`,
  },
  diatomics_run_time_sec: {
    key: `diatomics_run_time_sec`,
    property: `run_time_sec`,
    label: `Speed`,
    description: `Wall time in seconds for the full homonuclear diatomic curve sweep (H-U, 119 separations each), including calculator setup; summed over shards for parallel runs. Hardware varies by submission (recorded in the model YAML's hardware field); blank for submissions without recorded timings`,
    unit: `s`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~s`,
  },
  diatomics_time_multiplier: {
    key: `diatomics_time_multiplier`,
    label: `Slowdown`,
    description: `Diatomics wall time as a multiple of the fastest finite wall time among models matching the active task, energy-only, training-data and openness filters (1× = fastest model in view)`,
    unit: `×`,
    better: `lower`,
    format: `.2~f`,
  },
  // memory columns hidden by default until submissions record them (see MD equivalents)
  diatomics_max_gpu_mem_gb: {
    key: `diatomics_max_gpu_mem_gb`,
    property: `max_gpu_mem_gb`,
    label: `VRAM`,
    description: `Peak GPU memory (torch CUDA allocator high-water mark) over the full homonuclear diatomic sweep; max over shards for parallel runs`,
    unit: `GB`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~s`,
    visible: false,
  },
  diatomics_max_rss_gb: {
    key: `diatomics_max_rss_gb`,
    property: `max_rss_gb`,
    label: `RAM`,
    description: `Peak host memory (resident set size high-water mark) over the full homonuclear diatomic sweep; max over shards for parallel runs`,
    unit: `GB`,
    path: `metrics.diatomics`,
    better: `lower`,
    format: `.3~s`,
    visible: false,
  },
}

// Phonon metrics
export const PHONON_METRICS = {
  κ_SRME: {
    key: `κ_SRME`,
    label: `κ<sub>SRME</sub>`,
    description: `Symmetric relative mean error in predicted phonon-mode contributions to thermal conductivity κ, averaged over the 103 PhononDB-PBE materials (range [0, 2])`,
    path: `metrics.phonons.kappa_103`,
    range: [0, 2],
    better: `lower`,
    format: `.3~f`,
  },
  κ_SRE: {
    key: `κ_SRE`,
    label: `κ<sub>SRE</sub>`,
    description: `Symmetric relative error of total lattice thermal conductivity κ, averaged over the 103 PhononDB-PBE materials (range [0, 2])`,
    path: `metrics.phonons.kappa_103`,
    range: [0, 2],
    better: `lower`,
    format: `.3~f`,
  },
  κ_SRD: {
    key: `κ_SRD`,
    label: `κ<sub>SRD</sub>`,
    description: `Mean signed symmetric relative difference in total lattice thermal conductivity κ over the 103 PhononDB-PBE materials (range [-2, 2]); negative values indicate underprediction and positive values overprediction. Invalid or missing κ predictions count as zero conductivity (SRD = -2)`,
    path: `metrics.phonons.kappa_103`,
    range: [-2, 2],
    color_scale: `interpolateRdBu`,
    format: `.3~f`,
  },
  κ_failure_rate: {
    key: `κ_failure_rate`,
    label: `κ failed`,
    description: `Fraction of the 103 PhononDB-PBE materials where κ prediction failed outright and κ<sub>SRME</sub> was censored to its maximum of 2`,
    path: `metrics.phonons.kappa_103`,
    range: [0, 1],
    better: `lower`,
    format: `.1~%`,
  },
  imaginary_mode_rate: {
    key: `imaginary_mode_rate`,
    label: `Im(ω)`,
    description: `Fraction of the 103 PhononDB-PBE materials with imaginary phonon modes after ML relaxation; missing flags count as unflagged`,
    path: `metrics.phonons.kappa_103`,
    range: [0, 1],
    better: `lower`,
    format: `.1~%`,
  },
  spectrum_w1: {
    key: `spectrum_w1`,
    label: `W<sub>1</sub>(ω)`,
    description: `Mean Wasserstein-1 distance between ML and DFT phonon frequency spectra over materials with usable frequencies; robust to error compounding in κ, with an approximately 0.02 THz mesh-comparison noise floor`,
    path: `metrics.phonons.kappa_103`,
    unit: `THz`,
    better: `lower`,
    format: `.3~f`,
  },
} as const satisfies Record<string, Label>

type AllMetrics = DiscoveryMetricsLabels &
  GeoOptSymmetryMetricsLabels &
  MdMetricsLabels &
  Record<keyof typeof PHONON_METRICS, Label> &
  Record<DiatomicsMetricKey, Label> & {
    CPS: Label
    RMSD: Label
  }

export const ALL_METRICS: AllMetrics = {
  // Dynamic metrics
  CPS: {
    key: `CPS`,
    label: `CPS`,
    description: `Combined Performance Score averages discovery (F1), structure optimization (RMSD), and phonon performance (κ<sub>SRME</sub>) according to user-defined weights. Warning: This is not a stable metric. Further prediction tasks will be added to it in the future with the goal of making it a more holistic measure of overall model utility over time. When referring to it in papers, best include the benchmark version to avoid confusion (e.g. CPS v1 for the first version of CPS introduced in Matbench Discovery v1)`,
    range: [0, 1],
    better: `higher`,
    format: `.3f`,
  },
  ...DISCOVERY_METRICS,
  ...PHONON_METRICS,
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

// Column-visibility map for task pages: hide every metric except those passed in,
// keep all metadata columns visible. Keyed by col.key (unique per column) rather than
// label, which repeats across tasks (MD and diatomics both have Speed/Slowdown cols).
// Columns absent from the map (e.g. hyperparams) default to shown via `?? true` in
// the pages' col_filter.
export const task_page_visible_cols = (
  ...shown_metrics: Label[]
): Record<string, boolean> =>
  Object.fromEntries([
    ...Object.values(ALL_METRICS).map((col) => [col.key, false]),
    ...Object.values(METADATA_COLS).map((col) => [col.key, true]),
    ...shown_metrics.map((col) => [col.key, true]),
  ])

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

// Slowdown columns are roster-dependent (computed per filtered table view, not
// stored on models), so they can't be scatter axes
const time_multiplier_keys = new Set([
  DIATOMICS_METRICS.diatomics_time_multiplier.key,
  MD_METRICS.md_time_multiplier.key,
])
export const scatter_options = [
  ...Object.values(ALL_METRICS).filter((metric) => !time_multiplier_keys.has(metric.key)),
  HYPERPARAMS.model_params,
  METADATA_COLS.benchmark_added,
  METADATA_COLS.n_training_materials,
  METADATA_COLS.n_training_structures,
  HYPERPARAMS.graph_construction_radius,
  HYPERPARAMS.max_force,
  HYPERPARAMS.max_steps,
  HYPERPARAMS.batch_size,
  HYPERPARAMS.epochs,
  HYPERPARAMS.n_layers,
]

// Keyed lookup for bound axis/color selections.
export const scatter_options_by_key = Object.fromEntries(
  scatter_options.map((option) => [option.key, option]),
)

// Labels may contain HTML such as <sub>.
export const scatter_axis_label = (key: string): string =>
  scatter_options_by_key[key]?.label ?? key

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
// Add explicit mapping for hyperparams.
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
