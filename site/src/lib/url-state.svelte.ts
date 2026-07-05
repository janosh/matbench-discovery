import { afterNavigate, replaceState } from '$app/navigation'
import { page } from '$app/state'
import type { SortDir } from './types'

type PageState = Parameters<typeof replaceState>[1]
type SortState = { column: string; dir: SortDir }
type UrlParamEntry = [key: string, value: string, default_value?: string]
type ValidValues<T extends string> = ReadonlySet<T> | Record<string, unknown>
const sort_dirs = new Set<SortDir>([`asc`, `desc`])

export function valid_query_param<T extends string>(
  params: URLSearchParams,
  key: string,
  fallback: T,
  valid_values: ValidValues<T>,
): T {
  const value = params.get(key)
  if (!value) return fallback

  const is_valid =
    valid_values instanceof Set
      ? valid_values.has(value)
      : Object.hasOwn(valid_values, value)
  return is_valid ? (value as T) : fallback
}

// valid_columns is optional: the sortable-column set lives inside the table component
// (unknown columns are a harmless no-op there); pass it where known to reject garbage
export const sort_from_query = (
  params: URLSearchParams,
  default_sort: SortState,
  valid_columns?: ValidValues<string>,
): SortState => ({
  column: valid_columns
    ? valid_query_param(params, `sort`, default_sort.column, valid_columns)
    : (params.get(`sort`) ?? default_sort.column),
  dir: valid_query_param(params, `dir`, default_sort.dir, sort_dirs),
})

// -- Weighted-score (CPS/CMDS) radar weights as a single URL param --------------
// Serialized as comma-joined values in config-key order, e.g. weights=0.5,0.4,0.1.
type WeightsConfig = Record<string, { weight: number }>

const round_weight = (weight: number): number => Math.round(weight * 1000) / 1000

// Empty string when weights match the defaults (so sync_url_params drops the param)
export function weights_to_param(
  config: WeightsConfig,
  default_config: WeightsConfig,
): string {
  const keys = Object.keys(config)
  const is_default = keys.every(
    (key) =>
      round_weight(config[key].weight) ===
      round_weight(default_config[key]?.weight ?? NaN),
  )
  if (is_default) return ``
  return keys.map((key) => round_weight(config[key].weight)).join(`,`)
}

// Parse a weights param and write it into config (normalized to sum 1). A missing
// param resets to default_config: weight configs are shared module state, so without
// the reset, weights customized earlier in the session would survive navigating to a
// weights-less URL while all other URL-bound page state (sort, axes) resets.
// Silently ignores malformed params (wrong count, negative/non-finite, all-zero).
export function apply_weights_param(
  param: string | null,
  config: WeightsConfig,
  default_config: WeightsConfig,
): void {
  const keys = Object.keys(config)
  if (!param) {
    for (const key of keys) {
      config[key].weight = default_config[key]?.weight ?? config[key].weight
    }
    return
  }
  const values = param.split(`,`).map(Number)
  if (values.length !== keys.length) return
  if (values.some((val) => !Number.isFinite(val) || val < 0)) return
  const total = values.reduce((sum, val) => sum + val, 0)
  if (total === 0) return
  for (const [idx, key] of keys.entries()) config[key].weight = values[idx] / total
}

export function sync_url_params(entries: UrlParamEntry[], state: PageState): void {
  const params = new URLSearchParams(location.search)
  for (const [key, value, default_value = ``] of entries) {
    if (value === default_value) params.delete(key)
    else params.set(key, value)
  }

  const next_url = params.size ? `${location.pathname}?${params}` : location.pathname
  if (next_url !== `${location.pathname}${location.search}`) replaceState(next_url, state)
}

// Two-way URL query-param binding shared by all task pages. Reads state from the URL in
// afterNavigate (fires after the router is initialized, both on hydration and later
// navigations), then keeps the URL in sync with page state via replaceState. Gating
// writes on the first afterNavigate ensures the sync $effect never runs during the
// initial mount flush, which would throw "before router is initialized". Must be
// called during component init.
export function bind_url_params(
  read_params: ((params: URLSearchParams) => void) | null,
  entries: () => UrlParamEntry[],
): void {
  let url_ready = $state(false)

  afterNavigate(() => {
    read_params?.(page.url.searchParams)
    url_ready = true
  })

  $effect(() => {
    if (!url_ready) return
    sync_url_params(entries(), page.state)
  })
}
