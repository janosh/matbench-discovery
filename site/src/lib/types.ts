import type { Label as MattervizLabel } from 'matterviz'
import type { IconData } from 'svelte-widgets'
import type { Label1 as LabelType } from './schema/label'
import type { ModelMetadata, Person } from './schema/model'

export type { Dataset } from './schema/dataset'

export type OrgLogo =
  | { name: string; src: string; icon?: never }
  | { name: string; icon: IconData; src?: never }

export type ModelData = ModelMetadata & {
  // These fields are populated in MODELS variable in models.svelte.ts
  dirname: string
  metadata_file: string
  color?: string
  n_training_materials?: number
  n_training_structures?: number
  // training_cost summed over devices (count × hours per device), undefined if unreported
  training_gpu_hours?: number
  org_logos?: OrgLogo[]
  CPS?: number
}

export type Author = Person

export type Label = LabelType & Pick<MattervizLabel, `color_scale`>
export type TableLabel = Omit<Label, `better`> & MattervizLabel

export const DISCOVERY_SETS = [`full_test_set`, `unique_prototypes`] as const
export type DiscoverySet = (typeof DISCOVERY_SETS)[number]

export type SortDir = `asc` | `desc`

export interface DiatomicsCurves {
  distances: number[]
  'homo-nuclear': Record<string, { energies: number[]; forces: number[][] }>
  'hetero-nuclear'?: Record<string, { energies: number[]; forces: number[][] }>
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
