import { MD_METRICS } from '$lib/labels'
import MdPage from '$routes/tasks/md/+page.svelte'
import { describe, expect, it } from 'vitest'
import { doc_query, mount } from '../index'

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
    for (const header of [`ΔRDF`, `ΔADF`, `ΔvDOS`, `PMAE`, `PW1`, `CMDS`]) {
      expect(headers, `missing column ${header}`).toContain(header)
    }

    // DynamicScatter with h2 heading that tracks the selected axes
    const h2s = [...document.querySelectorAll<HTMLHeadingElement>(`h2`)]
    const scatter_heading = h2s.find(
      (h2) => h2.textContent?.replaceAll(/\s+/g, ` `).trim() === `ΔRDF vs FRMSE`,
    )
    expect(scatter_heading, `dynamic scatter heading`).toBeDefined()

    const scatter = doc_query<HTMLDivElement>(`div.scatter`)
    expect(scatter.getAttribute(`style`)).toContain(`height: 800px`)
  })

  it(`MD metric labels point at metrics.md; errors lower=better, CMDS higher=better`, () => {
    for (const label of Object.values(MD_METRICS)) {
      expect(label.path).toBe(`metrics.md`)
      // combined_score (CMDS) is a score (higher=better); the rest are errors (lower)
      const expected = label.key === `combined_score` ? `higher` : `lower`
      expect(label.better, label.key).toBe(expected)
    }
  })
})
