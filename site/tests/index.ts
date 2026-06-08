import { gzipSync } from 'node:zlib'
import { beforeAll, beforeEach, vi } from 'vitest'

// MatchMedia mock for Svelte MediaQuery - needed for svelte-multiselect
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

// Hoisted mocks for SvelteKit $app modules - more efficient than vi.mock at top level
const app_mocks = vi.hoisted(() => ({
  state: {
    page: {
      url: { pathname: `/`, searchParams: new URLSearchParams() },
      params: {},
      route: { id: null },
    },
  },
  environment: { browser: false, building: false, version: `test` },
  navigation: {
    goto: vi.fn(),
    invalidate: vi.fn(),
    invalidateAll: vi.fn(),
    preloadData: vi.fn(),
    preloadCode: vi.fn(),
    beforeNavigate: vi.fn(),
    afterNavigate: vi.fn(),
    pushState: vi.fn(),
    replaceState: vi.fn(),
  },
}))

vi.mock(`$app/state`, () => app_mocks.state)
vi.mock(`$app/environment`, () => app_mocks.environment)
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
})

// gzipped 200 Response for stubbing fetch() of .json.gz assets
export const gzipped_json_response = (data: unknown) =>
  Promise.resolve(new Response(gzipSync(JSON.stringify(data)), { status: 200 }))

// normalize the fetch() url argument (string | URL | Request) to a string
export const request_url = (url: RequestInfo | URL) =>
  typeof url === `string` || url instanceof URL ? String(url) : url.url

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
