import data from '$routes/models/per-element-each-errors.jsonl'

// Per-model mean convex-hull-distance error projected onto elements. Reference columns
// are payload-wide data in _base; only actual models occupy roster lines.
type ElementValues = Record<string, number | null>
type ElemErrModel = { model_key: string; values: ElementValues }
const { models, mp_occurrences, test_set_standard_deviation } = data as {
  models: ElemErrModel[]
  mp_occurrences: ElementValues
  test_set_standard_deviation: ElementValues
}
export const per_element_each_errors: Record<string, ElementValues> = Object.fromEntries([
  [`MP Occurrences`, mp_occurrences],
  [`Test set standard deviation`, test_set_standard_deviation],
  ...models.map((model) => [model.model_key, model.values] as const),
])
