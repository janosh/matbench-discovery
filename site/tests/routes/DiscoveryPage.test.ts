import { MODELS } from '$lib'
import { make_table_filters } from '$lib/models.svelte'
import DiscoveryPage from '$routes/tasks/discovery/+page.svelte'
import { tick } from 'svelte'
import { describe, expect, it } from 'vitest'
import {
  checkbox_for,
  filter_summary_badge,
  mount,
  mount_with_url,
  sorted_header,
} from '../index'

const active_toggle = (): string | undefined =>
  document.querySelector(`button.active`)?.textContent?.trim()

const button_for = (label: string): HTMLButtonElement => {
  const button = [...document.querySelectorAll<HTMLButtonElement>(`button`)].find(
    (candidate) => candidate.textContent?.trim() === label,
  )
  if (!button) throw new Error(`No button found for ${label}`)
  return button
}

const heading_texts = (): (string | undefined)[] =>
  [...document.querySelectorAll(`h2`)].map((heading) =>
    heading.textContent?.replaceAll(/\s+/g, ` `).trim(),
  )

describe(`Discovery Task Page`, () => {
  it(`renders page structure, discovery columns, and scatter`, () => {
    mount(DiscoveryPage, { target: document.body })

    expect(document.querySelector(`h1`)?.textContent).toContain(
      `Crystal Stability Prediction`,
    )

    const table = document.querySelector(`section.full-bleed table`)
    expect(table).not.toBeNull()

    const headers = [...document.querySelectorAll(`th .header-label`)].map((header) =>
      header.textContent?.trim(),
    )
    for (const col of [`Model`, `F1`, `DAF`, `Links`, `Date Added`]) {
      expect(headers).toContain(col)
    }
    expect(sorted_header()?.textContent).toContain(`F1`)
    expect(sorted_header()?.getAttribute(`aria-sort`)).toBe(`descending`)

    expect(heading_texts()).toContain(`F1 vs Params`)
    expect(document.querySelector(`[style*="height: 800px"]`)).not.toBeNull()

    expect(document.body.textContent).toContain(`Convex Hull Construction`)
  })

  it(`renders discovery-set controls and filters rows when toggled`, async () => {
    const default_filters = make_table_filters()
    const source_model = MODELS.find((model) => {
      const discovery = model.metrics?.discovery
      return (
        discovery != null &&
        typeof discovery === `object` &&
        default_filters.matches(model)
      )
    })
    const discovery = source_model?.metrics?.discovery
    if (!source_model || !discovery || typeof discovery !== `object`) {
      throw new Error(`No visible model with discovery metrics found`)
    }
    const partial_model = {
      ...source_model,
      model_key: `partial-discovery-test-model`,
      model_name: `Partial Discovery Test Model`,
      metrics: {
        ...source_model.metrics,
        discovery: { ...discovery, full_test_set: undefined },
      },
    }
    MODELS.push(partial_model)
    try {
      mount(DiscoveryPage, { target: document.body })
      await tick()

      const button_texts = [...document.querySelectorAll(`button`)].map((button) =>
        button.textContent?.trim(),
      )
      expect(button_texts).toEqual(
        expect.arrayContaining([`Unique Prototypes`, `Full Test Set`, `10k Most Stable`]),
      )
      expect(active_toggle()).toBe(`Unique Prototypes`)

      const table_text = () =>
        document.querySelector(`section.full-bleed tbody`)?.textContent
      expect(table_text()).toContain(partial_model.model_name)

      button_for(`Full Test Set`).click()
      await tick()

      expect(active_toggle()).toBe(`Full Test Set`)
      expect(table_text()).not.toContain(partial_model.model_name)
    } finally {
      MODELS.pop()
    }
  })

  it(`restores and syncs URL state`, async () => {
    const url = `http://localhost/tasks/discovery?set=full_test_set&targets=F,S,gradient&x=F1&y=rmsd&sort=F1&dir=asc&train=MPtrj,-OMat24&heatmap=0`
    await mount_with_url(DiscoveryPage, url)

    expect(active_toggle()).toBe(`Full Test Set`)
    expect(filter_summary_badge(`Targets`)).toContain(`(F,S,gradient)`)
    expect(heading_texts()).toContainEqual(expect.stringContaining(`RMSD vs F1`))
    expect(filter_summary_badge(`Training data`)).toContain(`(2)`)
    expect(checkbox_for(`Heatmap`).checked).toBe(false)
    const header = sorted_header()
    expect(header?.textContent).toContain(`F1`)
    expect(header?.getAttribute(`aria-sort`)).toBe(`ascending`)

    button_for(`10k Most Stable`).click()
    await tick()

    expect(new URL(location.href).searchParams.get(`set`)).toBe(`most_stable_10k`)
  })
})
