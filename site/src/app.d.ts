/// <reference types="@sveltejs/kit" />
/// <reference types="mdsvex/globals" />

declare module '*.md'

declare module '*package.json' {
  const pkg: Record<string, unknown>
  export default pkg
}

declare module '*metadata.yml' {
  const data: import('$lib').ModelMetadata
  export default data
}

declare module '*element-counts.json' {
  const map: Record<string, number>
  export default map
}
