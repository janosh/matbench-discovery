import {
  element_group_keys,
  element_groups,
  type ElementGroupKey,
} from 'matterviz/element'
import { DEFAULT_CATEGORY_COLORS } from 'matterviz/colors'
import type { SvelteSet } from 'svelte/reactivity'
import { in_view } from 'svelte-widgets/attachments'
import { valid_query_param, type UrlParamEntry } from 'svelte-widgets/url-params'

// Light, saturated-but-not-neon hues that read on both themes and carry black tooltip
// text. DFT references take cyan/periwinkle, which the model palette leaves out.
export const DFT_REF_COLORS: Readonly<Record<string, string>> = {
  PBE: `#22d3ee`,
  r2SCAN: `#818cf8`,
}

export const MODEL_COLORS = [
  `#fbbf24`, // amber
  `#34d399`, // emerald
  `#60a5fa`, // blue
  `#fb923c`, // orange
  `#a3e635`, // lime
  `#f87171`, // red
  `#2dd4bf`, // teal
  `#f9a8d4`, // pink
  `#d6d3d1`, // stone
] as const

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

// matterviz periodic-table category behind each element-group filter button (`all` has
// none). Nonmetals span two categories; halogens F-I are all diatomic nonmetals.
const GROUP_CATEGORIES: Partial<Record<ElementGroupKey, string>> = {
  alkali: `alkali metal`,
  alkaline_earth: `alkaline earth metal`,
  transition: `transition metal`,
  post_transition: `post-transition metal`,
  metalloid: `metalloid`,
  nonmetal: `polyatomic nonmetal`,
  halogen: `diatomic nonmetal`,
  noble_gas: `noble gas`,
  lanthanide: `lanthanide`,
  actinide: `actinide`,
}

// Attachment for an element-group ButtonGroup: tints each button with its category color
export const tint_element_groups = (node: HTMLElement): void => {
  for (const button of node.querySelectorAll<HTMLElement>(`button[data-value]`)) {
    const category = GROUP_CATEGORIES[button.dataset.value as ElementGroupKey]
    if (category) {
      button.style.setProperty(`--category-color`, DEFAULT_CATEGORY_COLORS[category])
    }
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
