export { default as AuthorBrief } from './AuthorBrief.svelte'
export { default as CaptionedMetricsTable } from './CaptionedMetricsTable.svelte'
export { default as Footer } from './Footer.svelte'
export { default as ModelCard } from './ModelCard.svelte'
export { default as Nav } from './Nav.svelte'
export { default as PtableHeatmap } from './PtableHeatmap.svelte'
export { default as PtableInset } from './PtableInset.svelte'
export { default as References } from './References.svelte'

export type ModelData = ModelMetadata & ModelStats

export type ModelMetadata = {
  model_name: string
  model_version: string
  matbench_discovery_version: string
  date_added: string
  date_published?: string
  authors: Author[]
  trained_by?: Author[]
  repo: string
  url?: string
  doi?: string
  paper?: string
  requirements?: Record<string, string>
  // whether this model was trained from scratch specifically for
  // this benchmark using our canonical training set
  trained_for_benchmark: boolean
  // whether model and/or data are open or closed
  open: 'OSOD' | 'CSOD' | 'OSCD' | 'CSCD'
  training_set: {
    title: string
    url: string
    n_structures: number
    n_materials?: number
  }
  hyperparams: Record<string, string | number>
  notes?: Record<string, string>
  dirname: string // models/{dirname}/{dirname}.yml
  n_params?: number // number of trainable parameters aka model size
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
  'Run Time (h)': number
  TPR: number
  TNR: number
  DAF: number
  GPUs: number
  CPUs: number
  slurm_jobs: number
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
  github?: string
}

// used in citation.cff
export type CffAuthor = Omit<Author, 'name'> & {
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
