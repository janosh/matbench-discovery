import { beforeAll, beforeEach, vi } from 'vitest'

beforeAll(() => {
  const animation = {
    pause: () => {},
    play: () => {},
    effect: {
      getComputedTiming: () => {
        return {}
      },
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
  parent: ParentNode = document,
): T {
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

// mock matchMedia browser API
globalThis.matchMedia = vi.fn()

// ResizeObserver mock
globalThis.ResizeObserver = class ResizeObserver {
  observe() {}
  unobserve() {}
  disconnect() {}
}

// TODO remove pending https://github.com/sveltejs/kit/issues/14143#issuecomment-3179138497
vi.stubGlobal(`__SVELTEKIT_PAYLOAD__`, { data: null })
