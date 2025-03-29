// This file is auto-generated from model-schema.yml. Do not edit directly.

export type ModelType = `GNN` | `UIP` | `BO-GNN` | `Fingerprint` | `Transformer` | `RF`
export type TargetType = `E` | `EF_G` | `EF_D` | `EFS_G` | `EFS_D` | `EFS_GM` | `EFS_DM`
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
  [k: string]: unknown
}
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

export interface ModelMetadata {
  model_name: string
  model_key?: string
  model_version: string
  matbench_discovery_version: string
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
  checkpoint_url: string | (`not available` | `missing`)
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
    radial_cutoff?: number
    [k: string]: unknown
  }
  notes?: {
    Description?: string
    Training?: string
    'Missing Preds'?: string
    html?: {
      [k: string]: unknown
    }
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
          [k: string]: unknown
        }
      }
    | `missing`
  train_task: `RP2RE` | `RS2RE` | `S2E` | `S2RE` | `S2EF` | `S2EFS` | `S2EFSM`
  test_task: `IP2E` | `IS2E` | `IS2RE` | `IS2RE-SR` | `IS2RE-BO`
  model_type: ModelType
  targets: TargetType
  openness: `OSOD` | `OSCD` | `CSOD` | `CSCD`
  status?: `aborted` | `complete` | `deprecated` | `pending` | `superseded`
  metrics?: {
    phonons?: PhononMetrics | (`not applicable` | `not available`)
    geo_opt?: GeoOptMetrics | (`not applicable` | `not available`)
    discovery?: DiscoveryMetrics
    diatomics?: DiatomicsMetrics | (`not applicable` | `not available`)
  }
}
export interface Person {
  name: string
  affiliation?: string
  email?: string
  url?: string
  orcid?: string
  github?: string
  corresponding?: boolean
}
export interface PhononMetrics {
  kappa_103?: {
    [k: string]: unknown
  }
  [k: string]: unknown
}
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
