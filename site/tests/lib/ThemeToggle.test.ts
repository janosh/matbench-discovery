import ThemeToggle from '$lib/ThemeToggle.svelte'
import { mount, tick } from 'svelte'
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest'

describe(`ThemeToggle`, () => {
  let prefers_dark: boolean

  beforeEach(() => {
    // Mock system preference
    prefers_dark = true
    Object.defineProperty(globalThis, `matchMedia`, {
      writable: true,
      value: vi.fn((query: string) => ({
        matches: prefers_dark,
        media: query,
        onchange: null,
        addListener: vi.fn(),
        removeListener: vi.fn(),
        addEventListener: vi.fn(),
        removeEventListener: vi.fn(),
        dispatchEvent: vi.fn(),
      })),
    })

    // Mock localStorage with property access support (e.g., localStorage.theme)
    const storage_map: Record<string, string> = {}
    const storage_base = {
      getItem: vi.fn((key: string) => storage_map[key] ?? null),
      setItem: vi.fn((key: string, value: string) => {
        storage_map[key] = value
      }),
      removeItem: vi.fn((key: string) => {
        delete storage_map[key]
      }),
      clear: vi.fn(() => {
        for (const key of Object.keys(storage_map)) delete storage_map[key]
      }),
    }

    const proxied_storage = new Proxy(
      storage_base as unknown as Storage & Record<string, string>,
      {
        get(target, prop, receiver) {
          if (typeof prop === `string` && prop in storage_map) return storage_map[prop]
          return Reflect.get(
            target as unknown as object,
            prop as unknown as PropertyKey,
            receiver as unknown as object,
          ) as unknown
        },
        set(target, prop, value, receiver) {
          if (typeof prop === `string`) {
            storage_map[prop] = String(value)
            return true
          }
          return Reflect.set(
            target as unknown as object,
            prop as unknown as PropertyKey,
            value as unknown,
            receiver as unknown as object,
          )
        },
      },
    )

    Object.defineProperty(globalThis, `localStorage`, {
      value: proxied_storage,
      writable: true,
    })
  })

  afterEach(vi.restoreAllMocks)

  it.each([
    {
      name: `uses saved theme from localStorage`,
      saved_theme: `light` as const,
      sys_dark: true,
      expected: `light`,
    },
    {
      name: `defaults to system preference (light)`,
      saved_theme: undefined,
      sys_dark: false,
      expected: `light`,
    },
    {
      name: `defaults to system preference (dark)`,
      saved_theme: undefined,
      sys_dark: true,
      expected: `dark`,
    },
  ])(`$name`, async ({ saved_theme, sys_dark, expected }) => {
    prefers_dark = sys_dark
    if (saved_theme) globalThis.localStorage.theme = saved_theme

    mount(ThemeToggle, { target: document.body })
    await tick()
    expect(document.documentElement.style.colorScheme).toBe(expected)
    expect(document.documentElement.dataset.theme).toBe(expected)
    const button_el = document.body.querySelector(`button`)
    const next = { light: `dark`, dark: `light` }[expected]
    expect(button_el?.getAttribute(`title`)).toContain(`Switch to ${next} theme`)
  })

  it(`toggles theme when clicked and updates colorScheme and data-theme`, () => {
    mount(ThemeToggle, { target: document.body })
    const button = document.body.querySelector(`button`)
    expect(document.documentElement.style.colorScheme).toBe(`dark`)
    expect(document.documentElement.dataset.theme).toBe(`dark`)
    expect(button?.getAttribute(`title`)).toContain(`Switch to light theme`)

    button?.click()
    expect(document.documentElement.style.colorScheme).toBe(`light`)
    expect(document.documentElement.dataset.theme).toBe(`light`)
    expect(globalThis.localStorage.theme).toBe(`light`)
    expect(button?.getAttribute(`title`)).toContain(`Switch to light theme`)

    button?.click()
    expect(document.documentElement.style.colorScheme).toBe(`dark`)
    expect(document.documentElement.dataset.theme).toBe(`dark`)
    expect(globalThis.localStorage.theme).toBe(`dark`)
  })
})
