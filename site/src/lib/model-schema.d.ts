// This file is auto-generated from model-schema.yml. Do not edit directly.

export type ModelMetadata = {
  [k: string]: unknown
} & {
  model_name: string
  model_key?: string
  model_version: string
  date_added: string
  date_published: string
  /**
   * @minItems 1
   */
  authors: {
    [k: string]: unknown
  } & [Person, ...Person[]]
  trained_by?: Person[]
  repo: string
  doi: string
  paper: string
  url?: string
  pypi?: string
  pr_url: string
  checkpoint_url: string | `missing`
  license: {
    /**
     * License for the model code
     */
    code:
      | `MIT`
      | `Apache-2.0`
      | `CC-BY-4.0`
      | `CC-BY-SA-4.0`
      | `CC-BY-NC-4.0`
      | `GPL-3.0`
      | `BSD-3-Clause`
      | `LGPL-3.0`
      | `Meta Research`
      | `ASL`
      | `unreleased`
    code_url?: string | `missing`
    /**
     * License for the model checkpoint
     */
    checkpoint:
      | `MIT`
      | `Apache-2.0`
      | `CC-BY-4.0`
      | `CC-BY-SA-4.0`
      | `CC-BY-NC-4.0`
      | `GPL-3.0`
      | `BSD-3-Clause`
      | `LGPL-3.0`
      | `Meta Research`
      | `ASL`
      | `unreleased`
    checkpoint_url?: string | `missing`
  }
  requirements: {
    /**
     * This interface was referenced by `undefined`'s JSON-Schema definition
     * via the `patternProperty` "^[a-zA-Z]{1}[a-zA-Z0-9_\-]{0,}$".
     */
    [k: string]: string
  }
  trained_for_benchmark: boolean
  training_set: (
    | `MP 2022`
    | `MPtrj`
    | `MPF`
    | `MP Graphs`
    | `GNoME`
    | `MatterSim`
    | `Alex`
    | `OMat24`
    | `sAlex`
    | `OpenLAM`
  )[]
  hyperparams?: {
    max_force?: number
    max_steps?: number
    optimizer?: string
    ase_optimizer?: string
    learning_rate?: number
    batch_size?: number
    epochs?: number
    n_layers?: number
    /**
     * Cutoff radius in Angstroms for graph construction (required for GNN and UIP models)
     */
    graph_construction_radius?: number
    /**
     * Maximum number of neighbors to consider in graph construction (required for GNN and UIP models)
     */
    max_neighbors?: number | `missing`
    [k: string]: unknown
  }
  notes?: {
    Description?: string
    Training?: string
    html?: string
    [k: string]: unknown
  }
  model_params: number
  n_estimators: number
  training_cost:
    | {
        /**
         * This interface was referenced by `undefined`'s JSON-Schema definition
         * via the `patternProperty` "^[a-zA-Z0-9\s]+ (GPUs|CPUs|TPUs)$".
         */
        [k: string]: {
          amount: number
          hours: number
          cost?: number
        }
      }
    | `missing`
  train_task: `RP2RE` | `RS2RE` | `S2E` | `S2RE` | `S2EF` | `S2EFS` | `S2EFSM`
  test_task: `IP2E` | `IS2E` | `IS2RE` | `IS2RE-SR`
  model_type: ModelType
  targets: TargetType
  openness: `OSOD` | `OSCD` | `CSOD` | `CSCD`
  status?: `aborted` | `complete` | `deprecated` | `pending` | `superseded`
  metrics?: {
    phonons?: PhononMetrics | (`not applicable` | `not available`)
    geo_opt?: GeoOptMetrics | (`not applicable` | `not available`)
    discovery?: DiscoveryMetrics | `not available`
    diatomics?: DiatomicsMetrics | (`not applicable` | `not available`)
  }
}
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "ModelType".
 */
export type ModelType = `GNN` | `UIP` | `BO-GNN` | `Fingerprint` | `Transformer` | `RF`
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "TargetType".
 */
export type TargetType = `E` | `EF_G` | `EF_D` | `EFS_G` | `EFS_D` | `EFS_GM` | `EFS_DM`
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "GeoOptMetrics".
 */
export type GeoOptMetrics = {
  [k: string]: unknown
} & {
  struct_col?: string
  pred_file?: string | null
  pred_file_url?: string
  'symprec=1e-5'?: {
    rmsd: number
    n_sym_ops_mae: number
    symmetry_decrease: number
    symmetry_match: number
    symmetry_increase: number
    n_structures: number
    analysis_file?: string | null
    analysis_file_url?: string
  }
  'symprec=1e-2'?: {
    rmsd: number
    n_sym_ops_mae: number
    symmetry_decrease: number
    symmetry_match: number
    symmetry_increase: number
    n_structures: number
    analysis_file?: string | null
    analysis_file_url?: string
  }
}
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "DiscoveryMetrics".
 */
export type DiscoveryMetrics = {
  [k: string]: unknown
} & {
  pred_col: string
  pred_file?: string | null
  pred_file_url?: string
  full_test_set?: DiscoveryMetricsSet
  most_stable_10k?: DiscoveryMetricsSet
  unique_prototypes?: DiscoveryMetricsSet
}
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "DiatomicsMetrics".
 */
export type DiatomicsMetrics = {
  [k: string]: unknown
} & {
  pred_file?: string | null
  pred_file_url?: string
  smoothness?: number
  tortuosity?: number
  conservation?: number
  energy_diff_flips?: number
  energy_grad_norm_max?: number
  energy_jump?: number
  force_mae?: number
  force_flips?: number
  force_total_variation?: number
  force_jump?: number
  force_conservation?: number
}
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "http_url".
 */
export type HttpUrl = string
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "pred_files".
 */
export type PredFiles = {
  [k: string]: unknown
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
  | `MIT`
  | `Apache-2.0`
  | `CC-BY-4.0`
  | `CC-BY-SA-4.0`
  | `CC-BY-NC-4.0`
  | `GPL-3.0`
  | `BSD-3-Clause`
  | `LGPL-3.0`
  | `Meta Research`
  | `ASL`
  | `unreleased`

/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "person".
 */
export interface Person {
  name: string
  affiliation?: string
  email?: string
  url?: string
  orcid?: string
  github?: string
  corresponding?: boolean
}
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "PhononMetrics".
 */
export interface PhononMetrics {
  kappa_103?: {
    [k: string]: unknown
  }
}
/**
 * This interface was referenced by `undefined`'s JSON-Schema
 * via the `definition` "DiscoveryMetricsSet".
 */
export interface DiscoveryMetricsSet {
  F1?: number
  DAF?: number
  Precision?: number
  Recall?: number
  Accuracy?: number
  TPR?: number
  FPR?: number
  TNR?: number
  FNR?: number
  TP?: number
  FP?: number
  TN?: number
  FN?: number
  MAE?: number
  RMSE?: number
  R2?: number
  missing_preds?: number
  missing_percent?: string
}
