/// <reference types="@sveltejs/kit" />
/// <reference types="mdsvex/globals" />

declare module '*.md'

declare module 'svelte-toc' {
  import type { Component, Snippet } from 'svelte'
  import type { HTMLAttributes } from 'svelte/elements'

  export interface TocProps extends HTMLAttributes<HTMLElement> {
    title?: string
    open?: boolean
    headings?: { id: string; text: string; level: number }[]
    headingSelector?: string
    sticky?: boolean
    activeHeading?: string
    activeHeadingScrollOffset?: number
    breakpoint?: number
    minItems?: number
    desktop?: boolean
    flashClickedHeadingsFor?: number
    getHeadingIds?: (node: HTMLElement) => string
    getHeadingLevels?: (node: HTMLElement) => number
    getHeadingText?: (node: HTMLElement) => string
    hide?: boolean
    keepActiveHeadingOnScroll?: boolean
    nav?: Snippet
    open_desktop?: Snippet
    open_mobile?: Snippet
    target_element?: HTMLElement
    warnOnEmpty?: boolean
    aside_style?: string
    nav_style?: string
  }

  const Toc: Component<TocProps>
  export default Toc
}

declare module 'svelte-toc/dist/MenuIcon.svelte' {
  import type { Component } from 'svelte'
  const MenuIcon: Component<Record<string, unknown>>
  export default MenuIcon
}

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
  export type DataFile = {
    url: string
    path: string
    description: string
    html?: string // auto-generated after ESM import in lib/index.ts
    figshare?: string
    md5?: string
  }
  type DataFiles = Omit<Record<string, DataFile>, '_links'> & {
    _links: string
  }
  const data_files: DataFiles
  export default data_files
}

declare module '*element-counts.json' {
  const map: Record<string, number>
  export default map
} // element counts for different datasets

declare module '*modeling-tasks.yml' {
  export type SubTask = {
    label: string
    description: string
  }

  export type ModelingTask = {
    label: string
    description: string
    metrics: {
      higher_is_better: string[]
      lower_is_better: string[]
    }
    subtasks?: Record<string, SubTask>
  }

  type ModelingTasks = {
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
