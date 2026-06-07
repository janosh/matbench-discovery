import data from '$routes/models/per-element-each-errors.json.gz'

// per-model mean convex-hull-distance error projected onto elements: column (a model_key
// or a metadata column like `Test set standard deviation`) -> element symbol -> value.
// gunzipped at build time by the json_gz plugin in vite.config.ts. Relative .json.gz
// imports only get the generic `unknown` type, so cast once here for all consumers.
export const per_element_each_errors = data as Record<
  string,
  Record<string, number | null>
>
