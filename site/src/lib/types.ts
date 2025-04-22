import * as d3sc from 'd3-scale-chromatic'
import { metric_labels } from './metrics'
import type { ModelMetadata } from './model-schema.d.ts'

export type { Dataset } from './dataset-schema.d.ts'
export type { ModelMetadata } from './model-schema.d.ts'

export type MetricKey = keyof typeof metric_labels
export type MetricLabel = (typeof metric_labels)[MetricKey]

export type ModelData = ModelMetadata &
  MetricKey & {
    // these fields are populated in MODELS variable in models.svelte.ts
    dirname: string
    metadata_file: string
    color?: string
    n_training_materials?: number
    n_training_structures?: number
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
  key: string // column header label
  label: string // column header label
  group?: string // group header label
  tooltip?: string // hover tooltip
  style?: string // CSS rules
  better?: `higher` | `lower` | null // sort direction
  color_scale?: keyof typeof d3sc // d3-scale-chromatic color scale name
  format?: string // d3-format string
  sticky?: boolean // sticky column
  visible?: boolean // show column (true by default)
  sortable?: boolean // whether column is sortable, defaults to true
  scale_type?: `linear` | `log` // scale type for color mapping
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
export type Metric = {
  path: string // path to the metric in the model metadata
  label: string // Display name (can include HTML)
  svg_label?: string // Label for the SVG chart
  description: string // Description of the metric
  range?: [number, number] // y-axis limits for the SVG chart
  better?: `higher` | `lower` // sort direction
}

export type CpsPart = `F1` | `kappa_SRME` | `RMSD`

export type CombinedMetricConfig = {
  parts: Record<CpsPart, Metric & { weight: number }>
  label: string
  name: string
  key: string
  range: [number, number]
  description: string
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
  | undefined
  | null
  | Record<string, unknown>
  | LinkData

export type RowData = {
  [key: string]: string | number | LinkData | null | undefined | boolean
}
