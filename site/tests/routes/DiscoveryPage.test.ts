import DiscoveryPage from '$routes/tasks/discovery/+page.svelte'
import { tick } from 'svelte'
import { describe, expect, it } from 'vitest'
import { mount, mount_with_url } from '../index'

const active_toggle = (): string | undefined =>
  document.querySelector(`button.active`)?.textContent?.trim()

const button_for = (label: string): HTMLButtonElement => {
  const button = [...document.querySelectorAll<HTMLButtonElement>(`button`)].find(
    (candidate) => candidate.textContent?.trim() === label,
  )
  if (!button) throw new Error(`No button found for ${label}`)
  return button
}

const energy_only_checkbox = (): HTMLInputElement => {
  const label = [...document.querySelectorAll(`label`)].find((candidate) =>
    candidate.textContent?.includes(`Energy-only models`),
  )
  const input = label?.querySelector<HTMLInputElement>(`input[type="checkbox"]`)
  if (!input) throw new Error(`Energy-only checkbox not found`)
  return input
}

const heading_texts = (): (string | undefined)[] =>
  [...document.querySelectorAll(`h2`)].map((heading) =>
    heading.textContent?.replaceAll(/\s+/g, ` `).trim(),
  )

describe(`Discovery Task Page`, () => {
  it(`renders page structure with title, table, and scatter`, () => {
    mount(DiscoveryPage, { target: document.body })

    expect(document.querySelector(`h1`)?.textContent).toContain(
      `Crystal Stability Prediction`,
    )

    const table = document.querySelector(`section.full-bleed table`)
    expect(table).not.toBeNull()

    expect(heading_texts()).toContain(`F1 vs Params`)
    expect(document.querySelector(`[style*="height: 800px"]`)).not.toBeNull()

    expect(document.body.textContent).toContain(`Convex Hull Construction`)
  })

  it(`renders SelectToggle with correct options and default selection`, () => {
    mount(DiscoveryPage, { target: document.body })

    const button_texts = [...document.querySelectorAll(`button`)].map((button) =>
      button.textContent?.trim(),
    )

    expect(button_texts).toContain(`Unique Prototypes`)
    expect(button_texts).toContain(`Full Test Set`)
    expect(button_texts).toContain(`10k Most Stable`)

    // Default selection
    expect(active_toggle()).toBe(`Unique Prototypes`)
  })

  it(`toggles discovery sets on button click`, async () => {
    mount(DiscoveryPage, { target: document.body })
    await tick()

    button_for(`Full Test Set`).click()
    await tick()

    expect(active_toggle()).toBe(`Full Test Set`)
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

  it(`restores and syncs URL state`, async () => {
    const url = `http://localhost/tasks/discovery?set=full_test_set&energy_only=1&x=F1&y=rmsd`
    await mount_with_url(DiscoveryPage, url)

    expect(active_toggle()).toBe(`Full Test Set`)
    expect(energy_only_checkbox().checked).toBe(true)
    expect(heading_texts()).toContainEqual(expect.stringContaining(`RMSD vs F1`))

    button_for(`10k Most Stable`).click()
    await tick()

    expect(new URL(location.href).searchParams.get(`set`)).toBe(`most_stable_10k`)
  })
})
