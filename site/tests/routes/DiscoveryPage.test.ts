import DiscoveryPage from '$routes/tasks/discovery/+page.svelte'
import { mount, tick } from 'svelte'
import { describe, expect, it } from 'vitest'

describe(`Discovery Task Page`, () => {
  it(`renders page structure with title, table, and scatter`, () => {
    mount(DiscoveryPage, { target: document.body })

    // Title
    expect(document.querySelector(`h1`)?.textContent).toContain(
      `Crystal Stability Prediction`,
    )

    // MetricsTable in full-bleed section
    const table = document.querySelector(`section.full-bleed table`)
    expect(table).not.toBeNull()

    // DynamicScatter with h2 heading that tracks the selected axes
    const h2s = [...document.querySelectorAll(`h2`)].map((h2) =>
      h2.textContent?.replaceAll(/\s+/g, ` `).trim(),
    )
    expect(h2s).toContain(`F1 vs Params`)
    expect(document.querySelector(`[style*="height: 800px"]`)).not.toBeNull()

    // Hull construction note from markdown
    expect(document.body.textContent).toContain(`Convex Hull Construction`)
  })

  it(`renders SelectToggle with correct options and default selection`, () => {
    mount(DiscoveryPage, { target: document.body })

    const buttons = document.querySelectorAll(`button`)
    const button_texts = [...buttons].map((btn) => btn.textContent?.trim())

    expect(button_texts).toContain(`Unique Prototypes`)
    expect(button_texts).toContain(`Full Test Set`)
    expect(button_texts).toContain(`10k Most Stable`)

    // Default selection
    expect(document.querySelector(`button.active`)?.textContent?.trim()).toBe(
      `Unique Prototypes`,
    )
  })

  it(`toggles discovery sets on button click`, async () => {
    mount(DiscoveryPage, { target: document.body })

    const full_test_btn = [...document.querySelectorAll(`button`)].find(
      (btn) => btn.textContent?.trim() === `Full Test Set`,
    )

    full_test_btn?.click()
    await tick()

    expect(document.querySelector(`button.active`)?.textContent?.trim()).toBe(
      `Full Test Set`,
    )
  })

  it(`shows discovery-specific columns in MetricsTable`, () => {
    mount(DiscoveryPage, { target: document.body })

    const headers = [...document.querySelectorAll(`th`)].map((th) =>
      th.textContent?.replace(/\s*[↑↓]\s*$/, ``).trim(),
    )

    for (const col of [`Model`, `F1`, `DAF`, `Links`, `Date Added`]) {
      expect(headers).toContain(col)
    }
  })
})
