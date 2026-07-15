import type * as d3sc from 'd3-scale-chromatic'
import type { IconName, Label as MattervizLabel } from 'matterviz'
import type { Label1 as LabelType } from './schema/label'
import type { ModelMetadata, Person } from './schema/model'

export type { Dataset } from './schema/dataset'

export type ModelData = ModelMetadata & {
  // These fields are populated in MODELS variable in models.svelte.ts
  dirname: string
  metadata_file: string
  color?: string
  n_training_materials?: number
  n_training_structures?: number
  org_logos?: { name: string; id?: string; src?: string; validated_icon?: IconName }[]
  CPS?: number
}

export type Author = Person

export type Label = LabelType & {
  color_scale?: keyof typeof d3sc // D3-scale-chromatic color scale name
  property?: string // Actual property name for data access (when different from key)
}
export type TableLabel = Omit<Label, `better`> & MattervizLabel

export const DISCOVERY_SETS = [
  `full_test_set`,
  `unique_prototypes`,
  `most_stable_10k`,
] as const
export type DiscoverySet = (typeof DISCOVERY_SETS)[number]

export type SortDir = `asc` | `desc`

export interface DiatomicsCurves {
  distances: number[]
  'homo-nuclear': Record<string, { energies: number[]; forces: number[][] }>
  'hetero-nuclear'?: Record<string, { energies: number[]; forces: number[][] }>
}

// Links data structure used for model resource links
export type LinkData = {
  paper: { url: string | null; title: string; icon: IconName }
  repo: { url: string | null; title: string; icon: IconName }
  pr_url: { url: string | null; title: string; icon: IconName }
  checkpoint?: { url: string | null; title: string; icon: IconName }
  pred_files: { files: { name: string; url: string }[]; name: string }
}

export interface GitHubActivityData {
  name: string
  repo: string
  stars: number
  forks: number
  commits_last_year: number
  contributors: number
  model_key?: string // URL slug for model detail page
}
