import type { ModelMetadata } from './model-schema.d.ts'

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
      n_structures: number
      n_materials?: number
      [k: string]: unknown
    }
