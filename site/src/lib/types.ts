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
}
