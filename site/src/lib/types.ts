import * as d3sc from 'd3-scale-chromatic'
import type { ModelMetadata } from './model-schema.d.ts'
export type * from './model-schema.d.ts'

export type ModelData = ModelMetadata &
  ModelStats & { dirname: string; metadata_file: string }
// dirname comes from: models/{dirname}/{model_name}.yml

export type ModelStats = {
  MAE: number // mean absolute error
  RMSE: number // root mean squared error
  R2: number
  Precision: number
  Recall: number
  F1: number
  missing_preds: number
  missing_percent: number
  Accuracy: number
  'Run Time (h)': number
  TPR: number // true positive rate
  TNR: number // true negative rate
  DAF: number // discovery acceleration factor
  GPUs: number // number of GPUs used
  CPUs: number // number of CPUs used
  slurm_jobs: number // number of SLURM jobs used
  κ_SRME: number // symmetric relative mean error for thermal conductivity
}

// how to pretty print a model stat key on the website
export type ModelStatLabel = {
  key: keyof ModelStats
  label?: string
  unit?: string
  tooltip?: string
}

export type Author = {
  name: string
  email?: string
  affiliation?: string
  orcid?: string
  url?: string
  twitter?: string
  github?: string
}

// used in citation.cff
export type CffAuthor = Omit<Author, `name`> & {
  'family-names': string
  'given-names': string
  affil_key: string
}

export type Reference = {
  title: string
  id: string
  author: { family: string; given: string }[]
  DOI: string
  URL?: string
  issued: { year: number; month: number; day: number }[]
  accessed: { year: number; month: number; day: number }[]
  page: string
  type: string
  ISSN?: string
}

export type Citation = {
  title: string
  subtitle?: string
  authors: CffAuthor[]
  affiliations: string[]
  'date-released': string
  license: string
  'license-url': string
  'repository-code': string
  url: string
  version: string
}

export type TrainingSet =
  | (`MP 2022` | `MPtrj` | `MPF` | `MP Graphs` | `GNoME` | `MatterSim` | `Alex`)
  | {
      title: string
      url: string
      download_url: string
      n_structures: number
      n_materials?: number
      [k: string]: unknown
    }

export type HeatmapColumn = {
  label: string // column header label
  group?: string // group header label
  tooltip?: string // hover tooltip
  style?: string // CSS rules
  better?: `higher` | `lower` | null // sort direction
  color_scale?: keyof typeof d3sc // d3-scale-chromatic color scale name
  format?: string // d3-format string
  sticky?: boolean // sticky column
  hidden?: boolean // hide column
  sortable?: boolean // whether column is sortable, defaults to true
}

export const DISCOVERY_SETS = [
  `full_test_set`,
  `unique_prototypes`,
  `most_stable_10k`,
] as const
export type DiscoverySet = (typeof DISCOVERY_SETS)[number]

export type DiatomicsCurves = {
  distances: number[]
  'homo-nuclear': Record<string, { energies: number[]; forces: number[][] }>
  'hetero-nuclear'?: Record<string, { energies: number[]; forces: number[][] }>
}

// MetricWeight defines weights for each component of the combined score
export interface MetricWeight {
  metric: string // ID of the metric (F1, RMSD, kappa_SRME)
  value: number // Weight value between 0-1
  label: string // Display name (can include HTML)
  description: string // Description of the metric
}

export interface CombinedMetricConfig {
  weights: MetricWeight[]
  name: string
  description: string
}

// Links data structure used for model resource links
export interface LinkData {
  paper: { url: string; title: string; icon: string }
  repo: { url: string; title: string; icon: string }
  pr_url: { url: string; title: string; icon: string }
  checkpoint?: { url: string | null; title: string; icon: string; is_missing?: boolean }
  pred_files: { files: { name: string; url: string }[]; name: string }
}

export type CellVal =
  | string
  | number
  | undefined
  | null
  | Record<string, unknown>
  | LinkData

export type RowData = Record<string, CellVal> & {
  row_attributes?: Record<string, string | boolean | number | undefined>
}

export type TableData = RowData[]
