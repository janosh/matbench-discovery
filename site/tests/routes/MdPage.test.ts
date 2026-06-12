import { MD_METRICS } from '$lib/labels'
import MdPage from '$routes/tasks/md/+page.svelte'
import { mount } from 'svelte'
import { describe, expect, it } from 'vitest'
import { doc_query } from '../index'

describe(`MD Task Page`, () => {
  it(`renders page structure with title, table, and scatter`, () => {
    mount(MdPage, { target: document.body })

    expect(document.querySelector(`h1`)?.textContent).toContain(
      `Molecular Dynamics Metrics`,
    )

    // MetricsTable in full-bleed section with MD metric columns visible
    const table = document.querySelector(`section.full-bleed table`)
    expect(table).not.toBeNull()
    const headers = [...document.querySelectorAll(`th`)].map((th) =>
      th.textContent?.replace(/\s*[↑↓]\s*$/, ``).trim(),
    )
    for (const header of [`RDF err`, `VDOS err`, `PMAE`, `PW1`, `Combined err`]) {
      expect(headers, `missing column ${header}`).toContain(header)
    }

    // DynamicScatter with h2 heading that tracks the selected axes
    const h2s = [...document.querySelectorAll<HTMLHeadingElement>(`h2`)]
    const scatter_heading = h2s.find(
      (h2) => h2.textContent?.replaceAll(/\s+/g, ` `).trim() === `RDF err vs FRMSE`,
    )
    expect(scatter_heading, `dynamic scatter heading`).toBeDefined()

    const scatter = doc_query<HTMLDivElement>(`div.scatter`)
    expect(scatter.getAttribute(`style`)).toContain(`height: 800px`)
  })

  it(`MD metric labels all point at metrics.md and are lower=better`, () => {
    for (const label of Object.values(MD_METRICS)) {
      expect(label.path).toBe(`metrics.md`)
      expect(label.better).toBe(`lower`)
    }
  })
})
