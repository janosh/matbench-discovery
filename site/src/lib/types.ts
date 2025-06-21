import * as d3sc from 'd3-scale-chromatic'
import type { Label1 as LabelType } from './label-schema.d.ts'
import type { AllMetrics } from './labels'
import type { ModelMetadata } from './model-schema.d.ts'

export type { Dataset } from './dataset-schema.d.ts'
export type { ModelMetadata } from './model-schema.d.ts'

export type ModelData = ModelMetadata &
  keyof AllMetrics & {
    // these fields are populated in MODELS variable in models.svelte.ts
    dirname: string
    metadata_file: string
    color?: string
    n_training_materials?: number
    n_training_structures?: number
    org_logos?: { name: string; id?: string; src?: string }[]
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

export type Label = LabelType & {
  color_scale?: keyof typeof d3sc // d3-scale-chromatic color scale name
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

// Links data structure used for model resource links
export type LinkData = {
  paper: { url: string; title: string; icon: string }
  repo: { url: string; title: string; icon: string }
  pr_url: { url: string; title: string; icon: string }
  checkpoint?: { url: string | null; title: string; icon: string; is_missing?: boolean }
  pred_files: { files: { name: string; url: string }[]; name: string }
}

export type CellVal =
  | string
  | number
  | boolean
  | undefined
  | null
  | Record<string, unknown>
  | LinkData
  | {
      [key: string]: string | number | LinkData | null | undefined | boolean
    }[]
export type RowData = { style?: string; [key: string]: CellVal }
export type CellSnippetArgs = { row: RowData; col: Label; val: CellVal }
