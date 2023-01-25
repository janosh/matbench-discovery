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
  run_time: number
  run_time_h: string
  GPUs: number
  CPUs: number
  slurm_jobs: number
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
