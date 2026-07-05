/// <reference types="@sveltejs/kit" />
/// <reference types="mdsvex/globals" />

declare module '*.md'

declare module '*.json.gz' {
  // data-only payloads, decompressed at build time by the json_payload plugin in
  // vite.config.ts. Figure payloads are typed in src/figs/payloads.d.ts.
  const data: unknown
  export default data
}

declare module '*.jsonl' {
  // line-delimited multi-model figure payloads, reassembled at build time by the
  // json_payload plugin in vite.config.ts, typed per payload in src/figs/payloads.d.ts
  const data: unknown
  export default data
}

declare module '*package.json' {
  const pkg: Record<string, unknown>
  export default pkg
}

declare module 'models/*.yml' {
  import type { ModelMetadata } from '$lib/schema/model'
  const data: ModelMetadata
  export default data
} // Model metadata files

declare module '*/datasets.yml' {
  import type { Dataset } from '$lib/types'
  const data: Record<string, Dataset>
  export default data
}

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
    md: ModelingTask
  }
  const tasks: ModelingTasks
  export default tasks
}

declare module '*mlip-github-activity.json' {
  import type { GitHubActivityData } from '$lib/types'
  const data: GitHubActivityData[]
  export default data
}
