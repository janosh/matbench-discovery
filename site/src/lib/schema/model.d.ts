// This file is auto-generated from model-schema.yml. Do not edit directly.
// oxlint-disable typescript/no-redundant-type-constituents

export type ModelMetadata = Record<string, unknown> & {
  model_name: string
  model_key: string
  /**
   * Upstream package / checkpoint version label (semver, date, tag, …).
   */
  model_version: string | null
  dates: Dates
  lifecycle: 'active' | 'superseded' | 'deprecated' | 'aborted'
  /**
   * model_key of the newer model replacing this one
   */
  superseded_by?: string
  /**
   * @minItems 1
   */
  authors: Record<string, unknown> & [Person, ...Person[]]
  trained_by?: Person[]
  repo: NullableHttpUrl
  doi: NullableHttpUrl
  paper: NullableHttpUrl
  /**
   * Model documentation / homepage URL (not the code repo).
   */
  docs?: HttpUrl | null
  pypi?: NullableHttpUrl
  pr_url: NullableHttpUrl
  checkpoint_url: NullableHttpUrl
  license: {
    /**
     * License for the model code
     */
    code:
      | 'MIT'
      | 'Apache-2.0'
      | 'CC-BY-4.0'
      | 'CC-BY-SA-4.0'
      | 'CC-BY-NC-4.0'
      | 'GPL-3.0'
      | 'BSD-3-Clause'
      | 'LGPL-3.0'
      | 'Meta Research'
      | 'ASL'
      | 'unreleased'
    code_url?: NullableHttpUrl
    code_url_reason?: string
    /**
     * License for the model checkpoint
     */
    checkpoint:
      | 'MIT'
      | 'Apache-2.0'
      | 'CC-BY-4.0'
      | 'CC-BY-SA-4.0'
      | 'CC-BY-NC-4.0'
      | 'GPL-3.0'
      | 'BSD-3-Clause'
      | 'LGPL-3.0'
      | 'Meta Research'
      | 'ASL'
      | 'unreleased'
    checkpoint_url?: NullableHttpUrl
    checkpoint_url_reason?: string
  }
  environment: Environment
  /**
   * Dataset keys validated against data/datasets.yml in registry tests.
   *
   * @minItems 1
   */
  training_sets: [string, ...string[]]
  /**
   * Optional model-size label (e.g. 10M, L, XL).
   */
  size?: string
  /**
   * Whether this checkpoint is a fine-tune of another training run.
   */
  fine_tune?: boolean
  hyperparams?: Hyperparams
  notes?: {
    description?: string
    training?: string
    html?: Record<string, string>
    [k: string]: unknown
  }
  model_params: number
  /**
   * Ensemble size. Omit for single models (default 1); only set when > 1.
   *
   */
  n_estimators?: number
  training_cost?: TrainingCost
  train_task: 'RP2RE' | 'RS2RE' | 'S2E' | 'S2RE' | 'S2EF' | 'S2EFS' | 'S2EFSM'
  test_task: 'IP2E' | 'IS2E' | 'IS2RE' | 'IS2RE-SR'
  /**
   * @minItems 1
   */
  architecture_types: [ArchitectureType, ...ArchitectureType[]]
  targets: TargetType
  openness: 'OSOD' | 'OSCD' | 'CSOD' | 'CSCD'
  metrics?: {
    phonons?: PhononMetrics
    geo_opt?: GeoOptMetrics
    discovery?: DiscoveryMetrics
    diatomics?: DiatomicsMetrics
    md?: MdMetrics
  }
}
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "http_url".
 */
export type HttpUrl = string
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "nullable_http_url".
 */
export type NullableHttpUrl = HttpUrl | null
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "ArchitectureType".
 */
export type ArchitectureType =
  | 'gnn'
  | 'transformer'
  | 'random_forest'
  | 'fingerprint'
  | 'bayesian_optimization'
  | 'unknown'
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "TargetType".
 */
export type TargetType =
  | 'E'
  | 'EF_G'
  | 'EF_D'
  | 'EFS_G'
  | 'EFSH_G'
  | 'EFS_D'
  | 'EFS_GM'
  | 'EFS_DM'
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "PhononMetrics".
 */
export type PhononMetrics = Record<string, unknown> & {
  status?: 'complete' | 'partial' | 'not_available' | 'not_applicable' | 'pending'
  reason?: string | null
  kappa_103?: KappaMetrics
}
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "nullable_file_ref".
 */
export type NullableFileRef = FileRef | null
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "GeoOptMetrics".
 */
export type GeoOptMetrics = Record<string, unknown> & {
  status?: 'complete' | 'partial' | 'not_available' | 'not_applicable' | 'pending'
  reason?: string | null
  pred_file?: NullableFileRef
  'symprec=1e-5'?: {
    rmsd: number
    n_sym_ops_mae: number
    symmetry_decrease: number
    symmetry_match: number
    symmetry_increase: number
    n_structures: number
    analysis_file?: NullableFileRef
  }
  'symprec=1e-2'?: {
    rmsd: number
    n_sym_ops_mae: number
    symmetry_decrease: number
    symmetry_match: number
    symmetry_increase: number
    n_structures: number
    analysis_file?: NullableFileRef
  }
}
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "DiscoveryMetrics".
 */
export type DiscoveryMetrics = Record<string, unknown> & {
  status?: 'complete' | 'partial' | 'not_available' | 'not_applicable' | 'pending'
  reason?: string | null
  pred_file?: NullableFileRef
  /**
   * Column name in pred_file containing formation-energy predictions. Prefer the canonical e_form_per_atom when present; required for older multi-column ensemble artifacts.
   *
   */
  pred_col?: string
  hardware?: string
  run_time_sec?: number
  max_rss_gb?: number
  max_gpu_mem_gb?: number
  full_test_set?: DiscoveryMetricsSet
  most_stable_10k?: DiscoveryMetricsSet
  unique_prototypes?: DiscoveryMetricsSet
}
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "DiatomicsMetrics".
 */
export type DiatomicsMetrics = Record<string, unknown> & {
  status?: 'complete' | 'partial' | 'not_available' | 'not_applicable' | 'pending'
  reason?: string | null
  pred_file?: NullableFileRef
  hardware?: string
  run_time_sec?: number
  max_rss_gb?: number
  max_gpu_mem_gb?: number
  excluded_formula_reasons?: Record<string, string>
  tortuosity?: number
  energy_diff_flips?: number
  energy_jump?: number
  pbe_wall_dist_mae?: number
  pbe_energy_mae?: number
  pbe_bond_length_error?: number
  pbe_well_depth_error?: number
  pbe_force_mae?: number
  pbe_vib_freq_error?: number
  force_flips?: number
  force_total_variation?: number
  force_jump?: number
}
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "MdMetrics".
 */
export type MdMetrics = Record<string, unknown> & {
  status?: 'complete' | 'partial' | 'not_available' | 'not_applicable' | 'pending'
  reason?: string | null
  pred_file?: NullableFileRef
  hardware?: string
  run_time_sec?: number
  max_rss_gb?: number
  max_gpu_mem_gb?: number
  energy_rmse?: number
  force_rmse?: number
  rdf_error?: number
  adf_error?: number
  vdos_error?: number
  pressure_mae?: number
  pressure_wasserstein?: number
  pressure_error?: number
  n_systems?: number
}
/**
 * License type:
 * - MIT: Massachusetts Institute of Technology
 * - CC-BY-4.0: Creative Commons Attribution 4.0 International
 * - CC-BY-NC-4.0: Creative Commons Attribution-NonCommercial 4.0 International
 * - CC-BY-SA-4.0: Creative Commons Attribution-ShareAlike 4.0 International
 * - GPL-3.0: GNU General Public v3.0
 * - BSD-3-Clause: Berkeley Software Distribution 3-Clause
 * - LGPL-3.0: GNU Lesser General Public License v3.0
 * - ASL: Academic Software License
 * - unreleased: No license since not released
 *
 *
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "license_enum".
 */
export type LicenseEnum =
  | 'MIT'
  | 'Apache-2.0'
  | 'CC-BY-4.0'
  | 'CC-BY-SA-4.0'
  | 'CC-BY-NC-4.0'
  | 'GPL-3.0'
  | 'BSD-3-Clause'
  | 'LGPL-3.0'
  | 'Meta Research'
  | 'ASL'
  | 'unreleased'

/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "Dates".
 */
export interface Dates {
  benchmark_added: string | null
  paper_published: string | null
}
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "person".
 */
export interface Person {
  name: string
  affiliation?: string
  email?: string
  url?: HttpUrl
  orcid?: HttpUrl
  github?: HttpUrl
  corresponding?: boolean
}
/**
 * Isolated runner / provenance environment. Keyed by model_key (dashes and dots mapped to underscores) when loading CALCULATORS specs.
 *
 *
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "Environment".
 */
export interface Environment {
  /**
   * Provenance dependency specs shown on the site. Without ``project``, these are also installed via ``uv run --with``. With ``project``, the project's pyproject.toml owns install resolution and these pins are display-only (keep them aligned with that project).
   *
   * @minItems 1
   */
  dependencies: [string, ...string[]]
  /**
   * Python version pin for the isolated uv runner (e.g. '3.12'). Use 3.11 only when a dependency stack requires it (e.g. older torch wheels).
   *
   */
  python_version: string
  /**
   * Extra ``uv --find-links`` URLs for the non-project runner path. Ignored when ``project`` is set (put find-links in that project's pyproject).
   *
   */
  find_links?: HttpUrl[]
  /**
   * Extra ``uv --extra-index-url`` entries for the non-project runner path. Ignored when ``project`` is set.
   *
   */
  extra_index_urls?: HttpUrl[]
  /**
   * Optional path to a uv project whose pyproject.toml owns dependency resolution (``uv run --project … --isolated``). Use when ``--with`` cannot express needed overrides (e.g. EquFlash).
   *
   */
  project?: string | null
}
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "Hyperparams".
 */
export interface Hyperparams {
  evaluation?: {
    /**
     * Geometry-optimization force convergence threshold in eV/Å.
     */
    max_force?: number
    force_max?: number
    max_steps?: number
    ase_optimizer?: string
    ase_filter?: string | null
    cell_filter?: string
    kappa?: KappaSettings
    [k: string]: unknown
  }
  architecture?: {
    /**
     * Graph construction cutoff radius in Å.
     */
    graph_construction_radius?: number
    /**
     * Maximum graph neighbors; null means unlimited.
     */
    max_neighbors?: number | null
    n_layers?: number
    [k: string]: unknown
  }
  training?: {
    learning_rate?: number
    initial_learning_rate?: number
    batch_size?: number
    epochs?: number
    /**
     * Canonical optimizer label (Adam, AdamW, Muon, ...).
     */
    optimizer?: string
    [k: string]: unknown
  }
  upstream_config?: Record<string, unknown>
}
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "KappaSettings".
 */
export interface KappaSettings {
  protocol: 'phonondb-v1'
  displacement_distance?: number
  /**
   * @minItems 1
   * @maxItems 1
   */
  temperatures?: [300]
  ase_optimizer?: string
  ase_filter?:
    | 'FrechetCellFilter'
    | 'ExpCellFilter'
    | 'UnitCellFilter'
    | 'frechet'
    | 'exp'
    | 'none'
    | null
  max_steps?: number
  force_max?: number
  symprec?: number
  relax_symprec?: number
  enforce_relax_symm?: boolean
  conductivity_broken_symm?: boolean
  is_plusminus?: true | false | 'auto'
  batch_size?: number
  max_atoms_per_batch?: number | null
  relaxation_mode?: 'single-stage' | 'two-stage' | 'none'
  save_forces?: boolean
}
/**
 * Reported training hardware/cost. Omit the whole field when unknown — do not add a placeholder reason.
 *
 *
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "TrainingCost".
 */
export interface TrainingCost {
  /**
   * @minItems 1
   */
  entries: [TrainingCostEntry, ...TrainingCostEntry[]]
  /**
   * Optional note about how cost was estimated or aggregated.
   */
  reason?: string
}
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "TrainingCostEntry".
 */
export interface TrainingCostEntry {
  /**
   * Accelerator / device label (e.g. NVIDIA H100).
   */
  hardware: string
  /**
   * Number of devices used.
   */
  count: number
  /**
   * Wall hours per device. Total device-hours = count × hours_per_device.
   */
  hours_per_device: number
}
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "KappaMetrics".
 */
export interface KappaMetrics {
  pred_file?: NullableFileRef
  force_file?: NullableFileRef
  run_info_file?: NullableFileRef
  hardware?: string
  run_time_sec?: number
  max_rss_gb?: number
  max_gpu_mem_gb?: number
  κ_SRME: number
  κ_SRE?: number
}
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "FileRef".
 */
export interface FileRef {
  name: string
  url?: HttpUrl
  size?: number
  md5?: string
}
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "DiscoveryMetricsSet".
 */
export interface DiscoveryMetricsSet {
  F1: number
  DAF: number
  Precision: number
  Recall: number
  Accuracy: number
  TPR: number
  FPR: number
  TNR: number
  FNR: number
  TP: number
  FP: number
  TN: number
  FN: number
  MAE: number
  RMSE: number
  R2: number
  missing_preds: number
}
