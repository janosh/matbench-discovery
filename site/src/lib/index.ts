export { default as Footer } from './Footer.svelte'
export { default as ModelCard } from './ModelCard.svelte'
export { default as Nav } from './Nav.svelte'
export { default as References } from './References.svelte'

export type ModelData = ModelMetadata & ModelStats

export type ModelMetadata = {
  model_name: string
  model_version: string
  matbench_discovery_version: string
  date_added: string
  authors: Author[]
  repo: string
  url?: string
  doi?: string
  preprint?: string
  requirements?: Record<string, string>
  trained_on_benchmark: boolean
}

export type ModelStats = {
  MAE: number
  RMSE: number
  R2: number
  Precision: number
  Recall: number
  F1: number
  missing_preds: number
  missing_percent: number
  Accuracy: number
  'Run Time (h)': string
  FPR: number
  FNR: number
  DAF: number
  GPUs: number
  CPUs: number
  slurm_jobs: number
  date_added: string
}

// [key, label?, unit?]
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
