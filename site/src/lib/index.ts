import type { ModelMetadata } from './model-schema'

export { default as CaptionedMetricsTable } from './CaptionedMetricsTable.svelte'
export { default as Footer } from './Footer.svelte'
export { default as AuthorBrief } from './ModelAuthor.svelte'
export { default as ModelCard } from './ModelCard.svelte'
export { default as Nav } from './Nav.svelte'
export { default as PtableHeatmap } from './PtableHeatmap.svelte'
export { default as PtableInset } from './PtableInset.svelte'
export { default as References } from './References.svelte'
export * from './types'

export function model_is_compliant(metadata: ModelMetadata): boolean {
  if ((metadata.openness ?? `OSOD`) != `OSOD`) return false

  const training_set = Array.isArray(metadata.training_set)
    ? metadata.training_set
    : [metadata.training_set]

  const allowed_sets = [`MP 2022`, `MPtrj`, `MPF`, `MP Graphs`]
  console.log(training_set)
  return training_set.every((itm) => allowed_sets.includes(itm))
}
