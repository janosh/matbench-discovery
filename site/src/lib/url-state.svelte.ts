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

export const sort_from_query = (
  params: URLSearchParams,
  default_sort: SortState,
): SortState => ({
  column: params.get(`sort`) ?? default_sort.column,
  dir: valid_query_param(params, `dir`, default_sort.dir, sort_dirs),
})

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
