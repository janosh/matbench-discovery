/// <reference types="@sveltejs/kit" />
/// <reference types="mdsvex/globals" />

declare module '*.md'
declare module '*package.json'

declare module '*metadata.yml' {
  const content: import('$lib/types').ModelMetadata
  export default content
}
