/// <reference types="@sveltejs/kit" />
/// <reference types="mdsvex/globals" />

declare module '*.md'

declare module '*package.json' {
  const pkg: Record<string, unknown>
  export default pkg
}

// model metadata files
declare module '*metadata.yml' {
  const data: import('$lib').ModelMetadata
  export default data
}

// paper metadata
declare module '*frontmatter.yml' {
  const frontmatter: import('$lib').Frontmatter
  export = frontmatter
}

declare module '*references.yaml' {
  export const references: import('$lib').Reference[]
}

declare module '*element-counts.json' {
  const map: Record<string, number>
  export default map
}
