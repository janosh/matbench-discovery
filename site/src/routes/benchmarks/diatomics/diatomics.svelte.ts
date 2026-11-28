import type { ChemicalElement, ElementCategory } from 'matterviz/element'
import { DEFAULT_CATEGORY_COLORS } from 'matterviz/colors'
import { SvelteMap, type SvelteSet } from 'svelte/reactivity'
import type { ButtonGroupOption } from 'svelte-widgets'
import { in_view } from 'svelte-widgets/attachments'
import { string_from_param, type UrlParamEntry } from 'svelte-widgets/browser'

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

// Element-group filter buttons. Elements belong by periodic-table `categories`, else by
// main-table `column` (halogens have no category of their own); `all` has neither. `tint` is the
// category whose color tints the button.
type ElementGroup = ButtonGroupOption & {
  categories?: readonly ElementCategory[]
  column?: number
  tint?: ElementCategory
}
const category_group = (value: string, label: string, category: ElementCategory) =>
  ({ value, label, categories: [category], tint: category }) satisfies ElementGroup

export const ELEMENT_GROUPS = [
  { value: `all`, label: `All`, tooltip: `Show all elements` },
  category_group(`alkali`, `Alkali metals`, `alkali metal`),
  category_group(`alkaline_earth`, `Alkaline earth metals`, `alkaline earth metal`),
  category_group(`transition`, `Transition metals`, `transition metal`),
  category_group(`post_transition`, `Post-transition metals`, `post-transition metal`),
  category_group(`metalloid`, `Metalloids`, `metalloid`),
  {
    value: `nonmetal`,
    label: `Nonmetals`,
    tooltip: `Diatomic and polyatomic nonmetals`,
    categories: [`diatomic nonmetal`, `polyatomic nonmetal`],
    tint: `polyatomic nonmetal`,
  },
  {
    value: `halogen`,
    label: `Halogens`,
    tooltip: `Group 17 halogen elements`,
    column: 17,
    tint: `diatomic nonmetal`, // F-I
  },
  category_group(`noble_gas`, `Noble gases`, `noble gas`),
  category_group(`lanthanide`, `Lanthanides`, `lanthanide`),
  category_group(`actinide`, `Actinides`, `actinide`),
] as const satisfies readonly ElementGroup[]

type ElementGroupKey = (typeof ELEMENT_GROUPS)[number][`value`]
const group_by_key = new Map<string, ElementGroup>(
  ELEMENT_GROUPS.map((group) => [group.value, group]),
)
const group_keys = new Set<ElementGroupKey>(ELEMENT_GROUPS.map(({ value }) => value))

// `elements` URL param: which element group's diatomics to show
export class ElementGroupParam {
  key = $state<ElementGroupKey>(`all`)
  readonly #group = $derived(group_by_key.get(this.key))

  includes = ({ category, column, row }: ChemicalElement): boolean => {
    if (!this.#group) throw new Error(`Unknown element group: ${this.key}`)
    const { categories, column: group_column } = this.#group
    if (categories) return categories.includes(category)
    // rows 8+ are the f-block, whose cells reuse main-table columns (Lu, Lr in 17)
    return group_column === undefined || (column === group_column && row <= 7)
  }

  read = (params: URLSearchParams): void => {
    this.key = string_from_param(params, `elements`, `all`, group_keys)
  }

  get url_entry(): UrlParamEntry {
    return [`elements`, this.key, `all`]
  }
}

// Attachment for an element-group ButtonGroup: tints each button with its category color
export const tint_element_groups = (node: HTMLElement): void => {
  for (const button of node.querySelectorAll<HTMLElement>(`button[data-value]`)) {
    const tint = group_by_key.get(button.dataset.value ?? ``)?.tint
    if (tint) button.style.setProperty(`--category-color`, DEFAULT_CATEGORY_COLORS[tint])
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

// === Curves served by curves/[key].json ===
// DFT references add per-formula distance grids plus the per-point electronic state
// shown in plot tooltips; model predictions carry energies only
export type PageDiatomicCurve = {
  energies: number[]
  distances?: number[]
  magmoms?: ([number, number] | null)[] // site-projected moments of the two atoms (μB)
  spin_candidates?: (string | null)[] // winning NUPDOWN (as string) or `afm`
  bond_forces?: number[] // force on atom 2 along the bond (eV/Å), > 0 = repulsive
}
export type PageDiatomicsCurves = {
  distances: number[]
  'homo-nuclear': Record<string, PageDiatomicCurve>
}

// Fetches each selected model's curves once; plots react as files arrive. A failed key
// is retried the next time it is requested.
export class DiatomicCurveCache {
  curves = new SvelteMap<string, PageDiatomicsCurves>()
  errors = new SvelteMap<string, string>()
  readonly #requested = new Set<string>()

  request = (keys: readonly string[]): void => {
    for (const key of keys) {
      if (this.#requested.has(key)) continue
      this.#requested.add(key)
      this.errors.delete(key)
      fetch(`/benchmarks/diatomics/curves/${encodeURIComponent(key)}.json`)
        .then(async (response) => {
          if (!response.ok) throw new Error(`HTTP ${response.status}`)
          this.curves.set(key, (await response.json()) as PageDiatomicsCurves)
        })
        .catch((error: unknown) => {
          this.#requested.delete(key)
          this.errors.set(key, error instanceof Error ? error.message : String(error))
        })
    }
  }
}
