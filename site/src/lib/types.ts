export type ModelMetadata = {
  model_name: string
  model_version: string
  matbench_discovery_version: string
  date_added: Date
  authors: Author[]
  repo: string
  url?: string
  doi?: string
  preprint?: string
  requirements?: Record<string, string>
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
