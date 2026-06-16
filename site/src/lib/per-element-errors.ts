import data from '$routes/models/per-element-each-errors.jsonl'

// per-model mean convex-hull-distance error projected onto elements: column (a model_key
// or a metadata column like `Test set standard deviation`) -> element symbol -> value.
// Stored line-delimited (one {key, values} column per line) so concurrent model
// submissions git-merge cleanly; reassembled by the jsonl Vite plugin, then flattened into
// the {column: {element: value}} map all consumers expect.
type ElemErrCol = { key: string; values: Record<string, number | null> }
const { models } = data as { models: ElemErrCol[] }
export const per_element_each_errors: Record<
  string,
  Record<string, number | null>
> = Object.fromEntries(models.map((col) => [col.key, col.values]))
