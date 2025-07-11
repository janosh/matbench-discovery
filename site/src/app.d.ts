/// <reference types="@sveltejs/kit" />
/// <reference types="mdsvex/globals" />

declare module '*.md'

declare module '*package.json' {
  const pkg: Record<string, unknown>
  export default pkg
}

declare module 'models/*.yml' {
  const data: import('$lib/model-schema').ModelMetadata
  export default data
} // model metadata files

declare module '*/datasets.yml' {
  const data: Record<string, import('$lib/types').Dataset>
  export default data
}

declare module '*citation.cff' {
  const data: import('$lib').Citation
  export default data
} // paper metadata

declare module '*references.yaml' {
  export const references: import('$lib').Reference[]
} // paper references (auto-exported by Zotero)

declare module '*model-schema.yml' {
  export const ModelMetadata: import('$lib/model-schema').ModelMetadata
} // model metadata schema

declare module '*data-files.yml' {
  type DataFile = {
    url: string
    path: string
    description: string
    html?: string // auto-generated after ESM import in lib/index.ts
    figshare?: string
    md5?: string
  }
  type DataFiles = {
    [K in string]: K extends `_links` ? string : DataFile
  }
  export const data_files: DataFiles
}

declare module '*element-counts.json' {
  const map: Record<string, number>
  export default map
} // element counts for different datasets

declare module '*modeling-tasks.yml' {
  type SubTask = {
    label: string
    description: string
  }

  type ModelingTask = {
    label: string
    description: string
    metrics: {
      higher_is_better: string[]
      lower_is_better: string[]
    }
    subtasks?: Record<string, SubTask>
  }
  export const geo_opt: ModelingTask
  export const discovery: ModelingTask
  export const phonons: ModelingTask
  export const diatomics: ModelingTask
}
