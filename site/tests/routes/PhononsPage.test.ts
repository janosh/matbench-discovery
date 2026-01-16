import PhononsPage from '$routes/tasks/phonons/+page.svelte'
import { mount } from 'svelte'
import { describe, expect, it } from 'vitest'

describe(`Phonons Task Page`, () => {
  it(`renders page structure with title, table, and scatter`, () => {
    mount(PhononsPage, { target: document.body })

    // Title
    expect(document.querySelector(`h1`)?.textContent).toContain(
      `MLFF Phonon Modeling Metrics`,
    )

    // MetricsTable in full-bleed section
    const table = document.querySelector(`section.full-bleed table`)
    expect(table).not.toBeNull()
    expect(document.querySelectorAll(`tbody tr`).length).toBeGreaterThan(0)

    // MetricScatter with h3 heading
    expect(document.querySelector(`h3`)?.textContent).toContain(`vs Model Parameters`)
    expect(document.querySelector(`[style*="height: 400px"]`)).not.toBeNull()
  })

  it(`shows κ_SRME and metadata columns, hides discovery metrics`, () => {
    mount(PhononsPage, { target: document.body })

    const headers = Array.from(document.querySelectorAll(`th`))
      .map((th) => th.textContent?.replace(/\s*[↑↓]\s*$/, ``).trim())

    // Should show κ_SRME and metadata
    expect(headers).toContain(`Model`)
    expect(headers).toContain(`Links`)
    expect(headers.some((header) => header?.includes(`κ`))).toBe(true)

    // Should hide discovery metrics
    for (const col of [`F1`, `DAF`, `Acc`, `Prec`]) {
      expect(headers).not.toContain(col)
    }
  })
})
