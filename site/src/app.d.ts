/// <reference types="@sveltejs/kit" />
/// <reference types="mdsvex/globals" />

declare module '*.md'

declare module '*package.json' {
  const pkg: Record<string, unknown>
  export default pkg
}

declare module 'models/*.yml' {
  const data: import('$lib/model-metadata').ModelMetadata
  export default data
} // model metadata files

declare module '*citation.cff' {
  const data: import('$lib').Citation
  export = data
} // paper metadata

declare module '*references.yaml' {
  export const references: import('$lib').Reference[]
} // paper references (auto-exported by Zotero)

declare module '*model-metadata-schema.yml' {
  export const ModelMetadata: import('$lib/model-metadata').ModelMetadata
} // model metadata schema

declare module '*element-counts.json' {
  const map: Record<string, number>
  export default map
}
