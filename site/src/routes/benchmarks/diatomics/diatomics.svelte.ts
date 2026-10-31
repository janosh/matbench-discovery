import {
  element_group_keys,
  element_groups,
  type ElementGroupKey,
} from 'matterviz/element'
import type { SvelteSet } from 'svelte/reactivity'
import { in_view } from 'svelte-widgets/attachments'
import { valid_query_param, type UrlParamEntry } from 'svelte-widgets/url-params'

// DFT reference colors outside the model palette, legible on light and dark themes
export const DFT_REF_COLORS: Readonly<Record<string, string>> = {
  PBE: `#17becf`,
  r2SCAN: `#f032e6`,
}

// `elements` URL param: which matterviz element group's diatomics to show
export class ElementGroupParam {
  key = $state<ElementGroupKey>(`all`)
  group = $derived(
    element_groups.find(({ value }) => value === this.key) ?? element_groups[0],
  )

  read = (params: URLSearchParams): void => {
    this.key = valid_query_param(params, `elements`, `all`, element_group_keys)
  }

  get url_entry(): UrlParamEntry {
    return [`elements`, this.key, `all`]
  }
}

// Attachment factory tracking which plot keys are near the viewport (render only those)
export const make_plot_observer =
  (visible: SvelteSet<string>) => (key: string) => (node: HTMLElement) => {
    const track = (is_visible: boolean) => visible[is_visible ? `add` : `delete`](key)
    const stop = in_view(track, { root_margin: `300px 0px` })(node)
    return () => {
      stop() // first, so a queued callback can't re-add the key
      visible.delete(key)
    }
  }
