/// <reference types="@sveltejs/kit" />
/// <reference types="mdsvex/globals" />

declare module '*.md'

declare module '*.json.gz' {
  // data-only figure payloads (site/src/figs/<name>.json.gz), decompressed at build
  // time by the json_gz_plugin in vite.config.ts and cast to shape in src/figs/index.ts
  const data: unknown
  export default data
}

declare module '*package.json' {
  const pkg: Record<string, unknown>
  export default pkg
}

declare module 'models/*.yml' {
  import type { ModelMetadata } from '$lib/model-schema'
  const data: ModelMetadata
  export default data
} // Model metadata files

declare module '*/datasets.yml' {
  import type { Dataset } from '$lib/types'
  const data: Record<string, Dataset>
  export default data
}

declare module '*citation.cff' {
  import type { Citation } from '$lib'
  const data: Citation
  export default data
} // Paper metadata

declare module '*references.yaml' {
  import type { Reference } from '$lib'
  export const references: Reference[]
} // Paper references (auto-exported by Zotero)

declare module '*model-schema.yml' {
  import type { ModelMetadata } from '$lib/model-schema'
  export const ModelMetadata: ModelMetadata
} // Model metadata schema

declare module '*data-files.yml' {
  export interface DataFile {
    url: string
    path: string
    description: string
    html?: string // Auto-generated after ESM import in lib/index.ts
    figshare?: string
    md5?: string
  }
  // Index signature allows DataFile entries, _links is the only string metadata field
  // Code filters keys starting with '_' when iterating over file entries
  interface DataFiles {
    [K: string]: DataFile | string
    _links: string
  }
  const data_files: DataFiles
  export default data_files
}

declare module '*element-counts.json' {
  const map: Record<string, number>
  export default map
} // Element counts for different datasets

declare module '*modeling-tasks.yml' {
  export interface SubTask {
    label: string
    description: string
  }

  export interface ModelingTask {
    label: string
    description: string
    metrics: {
      higher_is_better: string[]
      lower_is_better: string[]
    }
    subtasks?: Record<string, SubTask>
  }

  interface ModelingTasks {
    [key: string]: ModelingTask
    geo_opt: ModelingTask
    discovery: ModelingTask
    phonons: ModelingTask
    diatomics: ModelingTask
  }
  const tasks: ModelingTasks
  export default tasks
}

declare module '*mlip-github-activity.json' {
  import type { GitHubActivityData } from '$lib/types'
  const data: GitHubActivityData[]
  export default data
}
