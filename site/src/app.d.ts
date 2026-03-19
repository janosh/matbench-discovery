/// <reference types="@sveltejs/kit" />
/// <reference types="mdsvex/globals" />

declare module '*.md'

// Auto-generated Svelte figure components (excluded from TS checking in tsconfig.json)
declare module '$figs/*.svelte' {
  import type { Component } from 'svelte'
  const component: Component<{ style?: string; name?: string }>
  export default component
}

declare module '$figs/tmi/*.svelte' {
  import type { Component } from 'svelte'
  const component: Component<{ style?: string; name?: string }>
  export default component
}

declare module '*package.json' {
  const pkg: Record<string, unknown>
  export default pkg
}

declare module 'models/*.yml' {
  const data: import('$lib/model-schema').ModelMetadata
  export default data
} // Model metadata files

declare module '*/datasets.yml' {
  const data: Record<string, import('$lib/types').Dataset>
  export default data
}

declare module '*citation.cff' {
  const data: import('$lib').Citation
  export default data
} // Paper metadata

declare module '*references.yaml' {
  export const references: import('$lib').Reference[]
} // Paper references (auto-exported by Zotero)

declare module '*model-schema.yml' {
  export const ModelMetadata: import('$lib/model-schema').ModelMetadata
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
  const data: import('$lib/types').GitHubActivityData[]
  export default data
}
