import { beforeAll, beforeEach, vi } from 'vitest'

// matchMedia mock for Svelte MediaQuery - needed for svelte-multiselect
Object.defineProperty(window, `matchMedia`, {
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

// Mock $app modules for SvelteKit - must be at top level before any imports
vi.mock(`$app/state`, () => ({
  page: {
    url: {
      pathname: `/`,
      searchParams: new URLSearchParams(),
    },
    params: {},
    route: {
      id: null,
    },
  },
}))

vi.mock(`$app/environment`, () => ({ browser: false, building: false, version: `test` }))

vi.mock(`$app/navigation`, () => ({
  goto: vi.fn(),
  invalidate: vi.fn(),
  invalidateAll: vi.fn(),
  preloadData: vi.fn(),
  preloadCode: vi.fn(),
  beforeNavigate: vi.fn(),
  afterNavigate: vi.fn(),
  pushState: vi.fn(),
  replaceState: vi.fn(),
}))

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
