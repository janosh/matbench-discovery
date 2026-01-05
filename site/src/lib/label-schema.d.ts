// This file is auto-generated from label-schema.yml. Do not edit directly.

/**
 * This interface was referenced by `Label`'s JSON-Schema
 * via the `definition` "http_url".
 */
export type HttpUrl = string
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
 * This interface was referenced by `Label`'s JSON-Schema
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
 * This interface was referenced by `Label`'s JSON-Schema
 * via the `definition` "ModelType".
 */
export type ModelType = 'GNN' | 'UIP' | 'BO-GNN' | 'Fingerprint' | 'Transformer' | 'RF'
/**
 * This interface was referenced by `Label`'s JSON-Schema
 * via the `definition` "TargetType".
 */
export type TargetType = 'E' | 'EF_G' | 'EF_D' | 'EFS_G' | 'EFS_D' | 'EFS_GM' | 'EFS_DM'

export interface Label {
  [k: string]: unknown
}
/**
 * This interface was referenced by `Label`'s JSON-Schema
 * via the `definition` "Label".
 */
export interface Label1 {
  /**
   * Unique identifier for the metric
   */
  key: string
  /**
   * Short label for the metric
   */
  short?: string
  /**
   * Full label for the metric
   */
  label: string
  /**
   * Group header label
   */
  group?: string
  /**
   * Description of the metric
   */
  description: string
  /**
   * Path to the metric in the model metadata
   */
  path?: string
  /**
   * Unit of the metric
   */
  unit?: string
  /**
   * Range of possible values for the metric
   *
   * @minItems 2
   * @maxItems 2
   */
  range?: [number, number]
  /**
   * Whether the metric is visible
   */
  visible?: boolean
  /**
   * d3-format string
   */
  format?: string
  /**
   * Whether the metric is sticky
   */
  sticky?: boolean
  /**
   * Whether the metric is sortable
   */
  sortable?: boolean
  /**
   * Sort direction: higher is better, lower is better, or null for no sorting
   */
  better?: 'higher' | 'lower' | null
  /**
   * CSS string for cell style
   */
  cell_style?: string
  /**
   * CSS string for cell style
   */
  style?: string
  /**
   * Symprec value for the metric
   */
  symprec?: string
  /**
   * Scale type for the metric
   */
  scale_type?: 'log' | 'linear'
}
/**
 * This interface was referenced by `Label`'s JSON-Schema
 * via the `definition` "DiscoveryMetricsLabels".
 */
export interface DiscoveryMetricsLabels {
  Accuracy: Label1
  F1: Label1
  DAF: Label1
  Precision: Label1
  Recall: Label1
  TNR: Label1
  TPR: Label1
  MAE: Label1
  R2: Label1
  RMSE: Label1
}
/**
 * This interface was referenced by `Label`'s JSON-Schema
 * via the `definition` "MetadataLabels".
 */
export interface MetadataLabels {
  model_name: Label1
  training_set: Label1
  targets: Label1
  date_added: Label1
  links: Label1
  r_cut: Label1
  n_training_materials: Label1
  n_training_structures: Label1
  checkpoint_license: Label1
  code_license: Label1
  missing_preds: Label1
  'Run Time (h)': Label1
  org: Label1
}
/**
 * This interface was referenced by `Label`'s JSON-Schema
 * via the `definition` "HyperparamLabels".
 */
export interface HyperparamLabels {
  model_params: Label1
  graph_construction_radius: Label1
  max_force: Label1
  max_steps: Label1
  ase_optimizer: Label1
  cell_filter: Label1
  batch_size: Label1
  epochs: Label1
  n_layers: Label1
  learning_rate: Label1
  max_neighbors: Label1
  n_estimators: Label1
}
/**
 * This interface was referenced by `Label`'s JSON-Schema
 * via the `definition` "GeoOptSymmetryMetricsLabels".
 */
export interface GeoOptSymmetryMetricsLabels {
  'symmetry_match_1e-2': Label1
  'symmetry_decrease_1e-2': Label1
  'symmetry_increase_1e-2': Label1
  'symmetry_match_1e-5': Label1
  'symmetry_decrease_1e-5': Label1
  'symmetry_increase_1e-5': Label1
}
/**
 * This interface was referenced by `Label`'s JSON-Schema
 * via the `definition` "DatasetMetadataLabels".
 */
export interface DatasetMetadataLabels {
  name: Label1
  structures: Label1
  materials: Label1
  created: Label1
  open: Label1
  static: Label1
  license: Label1
  method: Label1
  api: Label1
  links: Label1
}
