import { beforeAll, beforeEach, vi } from 'vitest'

// matchMedia mock for Svelte MediaQuery - needed for svelte-multiselect
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

// Mock matterviz/table for ToggleMenu component
vi.mock(`matterviz/table`, async (importOriginal) => {
  const actual = await importOriginal<Record<string, unknown>>().catch(() => ({}))
  return {
    ...actual,
    ToggleMenu: () => null, // stub component
  }
})

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

export function doc_query<T extends HTMLElement>(
  selector: string,
  parent: ParentNode | null = document,
): T {
  if (!parent) throw new Error(`No parent node provided`)
  const node = parent.querySelector(selector)
  if (!node) throw new Error(`No element found for selector: ${selector}`)
  return node as T
}

// Helper function to check if an element is hidden
export function is_hidden(el: Element | null): boolean {
  if (!el) return true
  const style = getComputedStyle(el as HTMLElement)
  return (style.display === `none` || style.visibility === `hidden` ||
    el.getAttribute(`aria-hidden`) === `true` || el.hasAttribute(`hidden`))
}

// ResizeObserver mock
globalThis.ResizeObserver = class ResizeObserver {
  observe() {}
  unobserve() {}
  disconnect() {}
}
