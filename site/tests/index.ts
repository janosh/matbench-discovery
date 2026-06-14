import { gzipSync } from 'node:zlib'
import { mount as svelte_mount, unmount } from 'svelte'
import { afterEach, beforeAll, beforeEach, vi } from 'vitest'

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

// Svelte's mount() returns a live instance whose effects/subscriptions/listeners keep
// running until unmount(); clearing document.body.innerHTML only detaches the DOM. Without
// this, mounts leak across the file, accumulating effects that slow later renders until
// synchronous tests hit the 5s default timeout on CI. Track instances and unmount below.
// Typed concretely via Parameters/ReturnType (svelte-check-rs won't infer rest-arg types
// from a generic contextual annotation), then re-exported with svelte's generic mount
// signature so call sites keep full prop type-checking.
const mounted_components: ReturnType<typeof svelte_mount>[] = []
const tracked_mount = (
  ...args: Parameters<typeof svelte_mount>
): ReturnType<typeof svelte_mount> => {
  const instance = svelte_mount(...args)
  mounted_components.push(instance)
  return instance
}
export const mount = tracked_mount as typeof svelte_mount

afterEach(async () => {
  const instances = mounted_components.splice(0)
  await Promise.all(instances.map((instance) => unmount(instance)))
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
