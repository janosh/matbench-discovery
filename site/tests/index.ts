import { gzipSync } from 'node:zlib'
import { clear_asset_cache } from '#lib/asset-loader.js'
import type { Reproducibility } from '#lib/schema/model.js'
import type { ModelData } from '#lib/types.js'
import type { AfterNavigate } from '$app/navigation'
import NodeUtility from 'happy-dom/lib/nodes/node/NodeUtility.js'
import { mount as svelte_mount, tick, unmount } from 'svelte'
import {
  afterEach,
  beforeAll,
  beforeEach,
  expect,
  onTestFinished,
  vi,
} from 'vite-plus/test'

type AfterNavigateCallback = (navigation: unknown) => void

let after_navigate_callbacks: AfterNavigateCallback[] = []

// MatchMedia mock for Svelte MediaQuery - needed for svelte-widgets
Object.defineProperty(globalThis, `matchMedia`, {
  writable: true,
  value: (query: string) => ({
    matches: false,
    media: query,
    onchange: null,
    addListener: () => {},
    removeListener: () => {},
    addEventListener: () => {},
    removeEventListener: () => {},
    dispatchEvent: () => true,
  }),
})

// happy-dom lacks the Popover API that svelte-widgets >=1.6 tooltips use for their
// top-layer surface. Minimal polyfill: open state lives in a `data-popover-open` attribute
// (happy-dom can't match `:popover-open`, so `matches()` rewrites it to that attribute),
// and show/hide dispatch `toggle` so anything listening for popover state still hears it.
export const POPOVER_OPEN_ATTR = `data-popover-open`
if (!(`showPopover` in HTMLElement.prototype)) {
  const set_popover_open = (element: HTMLElement, open: boolean) => {
    if (!element.hasAttribute(`popover`)) {
      throw new DOMException(`Element is not a popover`, `NotSupportedError`)
    }
    if (element.hasAttribute(POPOVER_OPEN_ATTR) === open) return
    element.toggleAttribute(POPOVER_OPEN_ATTR, open)
    element.dispatchEvent(new Event(`toggle`))
  }
  Object.assign(HTMLElement.prototype, {
    showPopover(this: HTMLElement) {
      set_popover_open(this, true)
    },
    hidePopover(this: HTMLElement) {
      set_popover_open(this, false)
    },
    togglePopover(this: HTMLElement, force?: boolean) {
      const open = force ?? !this.hasAttribute(POPOVER_OPEN_ATTR)
      set_popover_open(this, open)
      return open
    },
  })
  const native_matches = Element.prototype.matches
  Object.defineProperty(Element.prototype, `matches`, {
    configurable: true,
    writable: true,
    value(this: Element, selectors: string) {
      return native_matches.call(
        this,
        selectors.replaceAll(`:popover-open`, `[${POPOVER_OPEN_ATTR}]`),
      )
    },
  })
}

// happy-dom can't match `:checked` on <option>, which Svelte's <select bind:value> queries
// on change, so resolve it from the select's own selection state
// oxlint-disable-next-line typescript/no-deprecated -- flags the deprecated-tag overloads, not the methods
const { querySelector, querySelectorAll } = Element.prototype
Object.assign(HTMLSelectElement.prototype, {
  querySelector(this: HTMLSelectElement, selectors: string) {
    if (selectors === `:checked`) return this.selectedOptions[0] ?? null
    return querySelector.call(this, selectors)
  },
  querySelectorAll(this: HTMLSelectElement, selectors: string) {
    if (selectors === `:checked`) return [...this.selectedOptions]
    return querySelectorAll.call(this, selectors)
  },
})

// happy-dom's IntersectionObserver never reports: treat every observed target as on
// screen so viewport-gated content (e.g. per-element diatomic plots) renders
globalThis.IntersectionObserver = class IntersectionObserver {
  readonly root = null
  readonly rootMargin = `0px`
  readonly thresholds = [0]
  readonly scrollMargin = `0px`
  constructor(private readonly callback: IntersectionObserverCallback) {}
  observe(target: Element): void {
    this.callback([{ isIntersecting: true, target } as IntersectionObserverEntry], this)
  }
  unobserve(): void {}
  disconnect(): void {}
  takeRecords = (): IntersectionObserverEntry[] => []
}

// happy-dom stubs `nodeName` as '' on Node.prototype and overrides it per subclass
// (Element, Text, ...). DOMPurify >= 3.4 reads nodeName through the Node.prototype getter
// it cached at construction (clobber-safe), so under happy-dom every node looks nameless
// and gets stripped, e.g. matterviz's sanitize_html(`r<sub>cut</sub>`) -> `cut`. Route the
// base getter to the nearest subclass getter so table headers keep their leading text.
if (Object.getOwnPropertyDescriptor(Node.prototype, `nodeName`)?.get?.call({}) === ``) {
  Object.defineProperty(Node.prototype, `nodeName`, {
    configurable: true,
    get(this: Node): string {
      let proto: object | null = Object.getPrototypeOf(this)
      while (proto && proto !== Node.prototype) {
        const getter = Object.getOwnPropertyDescriptor(proto, `nodeName`)?.get
        if (getter) return getter.call(this)
        proto = Object.getPrototypeOf(proto)
      }
      return ``
    },
  })
}

// happy-dom's NodeUtility.isFollowing (behind every Range boundary comparison, e.g. in
// extractContents) walks all nodes after node_b in document order until it meets node_a,
// each nextSibling step an O(siblings) indexOf. So one Range op costs O(document size), and
// matterviz's HeatmapTable, which splits every long HTML cell with a Range, renders full
// tables in quadratic time (~1/3 of MetricsTable test CPU). Same tree-order answer by
// comparing ancestor chains (equivalence checked in tests/setup.test.ts).
type HappyDomNode = Parameters<typeof NodeUtility.isFollowing>[0]
const inclusive_ancestors = (node: HappyDomNode): HappyDomNode[] => {
  const chain: HappyDomNode[] = []
  for (let current: HappyDomNode | null = node; current; current = current.parentNode) {
    chain.push(current)
  }
  return chain.toReversed() // tree root first
}
NodeUtility.isFollowing = (node_a: HappyDomNode, node_b: HappyDomNode): boolean => {
  if (node_a === node_b) return false
  const chain_a = inclusive_ancestors(node_a)
  const chain_b = inclusive_ancestors(node_b)
  if (chain_a[0] !== chain_b[0]) return false // separate trees: the walk never reaches node_a
  let depth = 1
  while (chain_a[depth] && chain_a[depth] === chain_b[depth]) depth++
  if (depth === chain_a.length) return false // node_a is an ancestor of node_b
  if (depth === chain_b.length) return true // node_a is a descendant of node_b
  const siblings = Array.from(chain_a[depth - 1].childNodes)
  return siblings.indexOf(chain_a[depth]) > siblings.indexOf(chain_b[depth])
}

// Hoisted mocks for SvelteKit $app modules - more efficient than vi.mock at top level
const app_mocks = vi.hoisted(() => ({
  state: {
    page: {
      url: new URL(`http://localhost/`),
      state: {},
    },
  },
  env: { building: false },
  navigation: {
    // shallow goto() only rewrites the URL, mirroring Kit's history update
    goto: vi.fn(
      async (url: string | URL, opts?: { shallow?: boolean; state?: unknown }) => {
        if (opts?.shallow) history.replaceState(opts.state, ``, url)
      },
    ),
    onNavigate: vi.fn(),
    afterNavigate: vi.fn((callback: AfterNavigateCallback) => {
      after_navigate_callbacks.push(callback)
    }),
  },
}))

vi.mock(`$app/state`, () => app_mocks.state)
vi.mock(`$app/env`, () => app_mocks.env)
vi.mock(`$app/navigation`, () => app_mocks.navigation)

beforeAll(() => {
  const animation = {
    pause: () => {},
    play: () => {},
    effect: {
      getComputedTiming: () => ({}),
      getKeyframes: () => [],
    },
    cancel: () => {},
    currentTime: 0,
  } as unknown as Animation

  Element.prototype.animate = () => animation
  Element.prototype.getAnimations = () => [animation]
})

beforeEach(() => {
  document.body.innerHTML = ``
  clear_asset_cache()
  const fetch_mock = vi.fn(async () => Response.error())
  vi.stubGlobal(`fetch`, fetch_mock)
  app_mocks.state.page.url = new URL(`http://localhost/`)
  after_navigate_callbacks = []
  history.replaceState(null, ``, `/`)
})

// Svelte's mount() returns a live instance whose effects/subscriptions/listeners keep
// running until unmount(); clearing document.body.innerHTML only detaches the DOM. Without
// this, mounts leak across the file, accumulating effects that slow later renders until
// synchronous tests hit the 5s default timeout on CI. Track instances and unmount below.
// Re-export with Svelte's generic mount signature so call sites retain prop type-checking.
const mounted_components: ReturnType<typeof svelte_mount>[] = []
const tracked_mount = (
  ...args: Parameters<typeof svelte_mount>
): ReturnType<typeof svelte_mount> => {
  const callback_start_idx = after_navigate_callbacks.length
  const instance = svelte_mount(...args)
  mounted_components.push(instance)
  queueMicrotask(() => {
    for (const callback of after_navigate_callbacks.slice(callback_start_idx)) {
      callback({ type: `enter`, from: null, to: app_mocks.state.page.url })
    }
  })
  return instance
}
export const mount = tracked_mount as typeof svelte_mount

export async function mount_with_url(
  component: unknown,
  url: string,
  options: Record<string, unknown> = {},
): Promise<ReturnType<typeof svelte_mount>> {
  const next_url = new URL(url)
  app_mocks.state.page.url = next_url
  history.replaceState(null, ``, `${next_url.pathname}${next_url.search}${next_url.hash}`)
  const instance = tracked_mount(component as Parameters<typeof svelte_mount>[0], {
    ...options,
    target: document.body,
  })
  await tick()
  return instance
}

export const query_param = (key: string) => new URLSearchParams(location.search).get(key)

export async function navigate(
  url: string,
  type: AfterNavigate[`type`] = `link`,
): Promise<void> {
  const destination = new URL(url, location.origin)
  app_mocks.state.page.url = destination
  history.replaceState(
    null,
    ``,
    `${destination.pathname}${destination.search}${destination.hash}`,
  )
  for (const [callback] of app_mocks.navigation.afterNavigate.mock.calls)
    callback({ type, from: null, to: { url: destination } })
  await tick()
}

export const sorted_header = (): HTMLTableCellElement | null =>
  document.querySelector(`thead th[aria-sort]:not([aria-sort="none"])`)

// column label alone: HeatmapTable appends a sort indicator (↓/↑ plus rank <sup>) and,
// for date columns, a format menu, both inside the <th>
export function header_name(header: Element): string {
  const label = header.cloneNode(true) as Element
  for (const control of label.querySelectorAll(`.header-popover, .resize-handle`)) {
    control.remove()
  }
  return (label.textContent ?? ``).trim().replace(/\s*[↓↑]\d*$/u, ``)
}

// table-controls filter dropdown trigger by its label (e.g. `Training data (2)`)
export function filter_menu_trigger(menu_name: string): HTMLButtonElement {
  const trigger = [
    ...document.querySelectorAll<HTMLButtonElement>(`button.filter-menu-trigger`),
  ].find((button) => button.textContent?.includes(menu_name))
  if (!trigger) throw new Error(`No filter menu found for ${menu_name}`)
  return trigger
}
export const filter_summary_badge = (menu_name: string): string =>
  filter_menu_trigger(menu_name).textContent?.trim() ?? ``

// open a filter dropdown (its Popover renders the panel only while open), return the panel
export async function open_filter_menu(menu_name: string): Promise<HTMLElement> {
  const trigger = filter_menu_trigger(menu_name)
  if (trigger.getAttribute(`aria-expanded`) !== `true`) {
    trigger.click()
    await tick()
  }
  const panel_id = trigger.getAttribute(`aria-controls`)
  const panel = document.querySelector<HTMLElement>(`[id="${panel_id}"]`)
  if (!panel) throw new Error(`Filter menu ${menu_name} did not open`)
  return panel
}

// open ColumnPicker's column popover (rendered only while open), return the panel
export async function open_column_menu(): Promise<HTMLElement> {
  const trigger = doc_query<HTMLButtonElement>(`.column-picker > button`)
  if (trigger.getAttribute(`aria-expanded`) !== `true`) {
    trigger.click()
    await tick()
  }
  return doc_query(`.column-picker .column-menu`)
}
export const heatmap_toggle = async (): Promise<HTMLInputElement> =>
  doc_query(`input[aria-label="Toggle heatmap colors"]`, await open_column_menu())

// a task page's metrics table restored `sort`/`dir` and `heatmap=0` from its URL;
// returns the heatmap checkbox for follow-up toggling
export async function expect_table_url_state(
  sort_label: string,
  aria_sort: `ascending` | `descending`,
): Promise<HTMLInputElement> {
  const header = sorted_header()
  if (!header) throw new Error(`No sorted table header`)
  expect(header_name(header)).toBe(sort_label)
  expect(header.getAttribute(`aria-sort`)).toBe(aria_sort)
  const heatmap = await heatmap_toggle()
  expect(heatmap.checked).toBe(false)
  return heatmap
}

// first element matching selector whose trimmed text equals text exactly
export function button_for(text: string, selector = `button`): HTMLButtonElement {
  const button = [...document.querySelectorAll<HTMLButtonElement>(selector)].find(
    (candidate) => candidate.textContent?.trim() === text,
  )
  if (!button) throw new Error(`No ${selector} found with text ${text}`)
  return button
}

// whitespace-normalized text of every heading matching selector (default h2)
export const heading_texts = (selector = `h2`): (string | undefined)[] =>
  [...document.querySelectorAll(selector)].map((heading) =>
    heading.textContent?.replaceAll(/\s+/g, ` `).trim(),
  )

// happy-dom lays nothing out, so plots need explicit client sizes to render markers
export function mock_element_size(width: number, height: number): void {
  vi.spyOn(HTMLElement.prototype, `clientWidth`, `get`).mockReturnValue(width)
  vi.spyOn(HTMLElement.prototype, `clientHeight`, `get`).mockReturnValue(height)
}

// find a checkbox by the text of its wrapping <label> (e.g. table-control toggles)
export function checkbox_for(label_text: string): HTMLInputElement {
  const label = [...document.querySelectorAll(`label`)].find((lbl) =>
    lbl.textContent?.includes(label_text),
  )
  const input = label?.querySelector<HTMLInputElement>(`input[type="checkbox"]`)
  if (!input) throw new Error(`No checkbox found for label ${label_text}`)
  return input
}

// Clear the site-wide comparison selection (module state shared by every test in a file)
// and close its dialog. Lazy import: a static one would load #lib before test files
// register their vi.mock()s.
export async function reset_comparison(): Promise<void> {
  const { comparison } = await import(`#lib/model-comparison.svelte.js`)
  comparison.keys.clear()
  comparison.open = false
}

// Shared task-page skeleton: section headings in order, an always-expanded (never
// collapsible) methodology last and, for pages with a combined score, its weights
export function expect_task_page_layout(headings: string[], n_score_weights = 0): void {
  expect(heading_texts()).toEqual(headings)
  expect(document.querySelector(`#methodology ~ details`)).toBeNull()
  expect(document.querySelectorAll(`.score-weights input[type="number"]`)).toHaveLength(
    n_score_weights,
  )
  expect(document.querySelector(`details .score-weights`)).toBeNull()
}

// Override fields of a shared module-level model (e.g. an entry of MODELS) for the
// current test only; fields the patch added are deleted again afterwards.
export function patch_model<Model extends object>(
  model: Model,
  patch: Partial<Model>,
): Model {
  const original = { ...model }
  Object.assign(model, patch)
  onTestFinished(() => {
    for (const key of Object.keys(patch))
      if (!(key in original)) Reflect.deleteProperty(model, key)
    Object.assign(model, original)
  })
  return model
}

// Synthetic metrics.<task>.reproducibility blocks, one per status
const repro = (patch: Partial<Reproducibility>): Reproducibility => ({
  checked: `2026-10-20`,
  status: `reproduced`,
  tolerance: 0.001,
  runs: [{ label: `H200 rerun`, file: { name: `rerun-preds.csv.gz` } }],
  rel_changes: { force_rmse: 0.0002 },
  ...patch,
})
export const REPRO_CASES = {
  reproduced: repro({}),
  drifted: repro({ status: `drifted`, cause: `environment_change`, note: `torch 2.9` }),
  inconclusive: repro({ status: `inconclusive` }),
}

// Add a synthetic model to a shared roster (e.g. ACTIVE_MODELS) for the current test only
export function push_model(roster: ModelData[], model: ModelData): ModelData {
  roster.push(model)
  onTestFinished(() => {
    roster.splice(roster.indexOf(model), 1)
  })
  return model
}

// trimmed labels of the MultiSelect chips (selected options) inside `scope`
export const selected_chips = (scope: ParentNode | string = document): string[] =>
  [
    ...(typeof scope === `string` ? doc_query(scope) : scope).querySelectorAll(
      `ul[aria-label="selected options"] > li`,
    ),
  ].map((chip) => chip.textContent?.trim() ?? ``)

// leaderboard toggle named only by its visible label (`Show all` / `Show only N selected`)
export const selected_only_toggle = (): HTMLInputElement | null =>
  [...document.querySelectorAll(`label`)]
    .find((label) =>
      /^Show (?:all|only \d+ selected)$/.test(label.textContent?.trim() ?? ``),
    )
    ?.querySelector<HTMLInputElement>(`input[type="checkbox"]`) ?? null

afterEach(async () => {
  const instances = mounted_components.splice(0)
  await Promise.all(instances.map((instance) => unmount(instance)))
  vi.restoreAllMocks()
})

// gzipped 200 Response for stubbing fetch() of .json.gz assets
export const gzipped_json_response = async (data: unknown) =>
  new Response(gzipSync(JSON.stringify(data)), { status: 200 })

// normalize the fetch() url argument (string | URL | Request) to a string
export const request_url = (url: RequestInfo | URL) =>
  typeof url === `string` || url instanceof URL ? String(url) : url.url

// stub fetch() with gzipped JSON for exact URLs and 404 for anything else
export const stub_json_routes = (routes: Record<string, unknown>) => {
  const fetch_mock = vi.fn((url: RequestInfo | URL) => {
    const href = request_url(url)
    return Object.hasOwn(routes, href)
      ? gzipped_json_response(routes[href])
      : Promise.resolve(new Response(`missing`, { status: 404 }))
  })
  vi.stubGlobal(`fetch`, fetch_mock)
  return fetch_mock
}

export const has_md_metrics = (model: ModelData): boolean => model.metrics?.md != null

export function doc_query<T extends Element = HTMLElement>(
  selector: string,
  parent: ParentNode | null = document,
  element_type?: abstract new (...args: never[]) => T,
): T {
  if (!parent) throw new Error(`No parent node provided`)
  const node = parent.querySelector(selector)
  if (!node) throw new Error(`No element found for selector: ${selector}`)
  if (element_type && !(node instanceof element_type)) {
    throw new Error(`Element for selector ${selector} has unexpected type`)
  }
  return node as T
}

export async function choose_scatter_property(
  label: string,
  query: string,
): Promise<void> {
  const is_size = label === `Marker size`
  const picker = is_size ? doc_query(`.property-picker`) : document
  if (is_size) {
    const input = doc_query<HTMLInputElement>(`input[role="combobox"]`, picker)
    input.focus()
    input.value = query
    input.dispatchEvent(new InputEvent(`input`, { bubbles: true }))
  } else {
    doc_query<HTMLButtonElement>(
      label === `Color`
        ? `.colorbar .property-select`
        : `.${label[0].toLowerCase()}-label button`,
    ).click()
  }
  await tick()
  const options = picker.querySelectorAll<HTMLElement>(
    is_size ? `ul.options li[aria-posinset]` : `.popover-select-dropdown [role="option"]`,
  )
  // size options render the (HTML) label followed by a <small> data-path/count line
  const option_text = (element: HTMLElement): string | undefined =>
    is_size
      ? element
          .querySelector(`span`)
          ?.textContent?.replace(element.querySelector(`small`)?.textContent ?? ``, ``)
      : element.textContent
  const option = [...options].find((element) => option_text(element)?.trim() === query)
  if (!option) throw new Error(`Missing option ${query} in ${label}`)
  option.click()
  await tick()
}

export function get_scatter_plot_props(scatter_plot_mock: {
  mock: { lastCall?: unknown[] }
}): unknown {
  const props = scatter_plot_mock.mock.lastCall?.find(
    (call_arg) =>
      typeof call_arg === `object` && call_arg !== null && `series` in call_arg,
  )
  if (!props) throw new Error(`ScatterPlot props not found`)
  return props
}

// ResizeObserver mock
globalThis.ResizeObserver = class ResizeObserver {
  observe() {}
  unobserve() {}
  disconnect() {}
}
