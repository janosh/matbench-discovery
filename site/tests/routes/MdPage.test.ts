import { MODELS } from '$lib'
import { MD_METRICS } from '$lib/labels'
import MdPage from '$routes/tasks/md/+page.svelte'
import { describe, expect, it } from 'vitest'
import { doc_query, has_md_metrics, mount, mount_with_url, sorted_header } from '../index'

describe(`MD Task Page`, () => {
  it(`renders page structure with filtered leaderboard and scatter`, () => {
    mount(MdPage, { target: document.body })

    expect(document.querySelector(`h1`)?.textContent).toContain(
      `Molecular Dynamics Metrics`,
    )

    const table = doc_query(`section.full-bleed table`)
    expect(table.querySelectorAll(`tbody tr`)).toHaveLength(
      MODELS.filter(has_md_metrics).length,
    )
    const headers = [...table.querySelectorAll(`th`)].map((th) =>
      th.textContent?.replace(/\s*[↑↓]\s*$/, ``).trim(),
    )
    for (const header of [`ΔRDF`, `ΔADF`, `ΔvDOS`, `PMAE`, `PW1`, `CMDS`]) {
      expect(headers, `missing column ${header}`).toContain(header)
    }

    const headings = [...document.querySelectorAll<HTMLHeadingElement>(`h2`)].map((h2) =>
      h2.textContent?.replaceAll(/\s+/g, ` `).trim(),
    )
    expect(headings).toContain(`ΔRDF vs FRMSE`)

    const scatter = doc_query<HTMLDivElement>(`div.scatter`)
    expect(scatter.getAttribute(`style`)).toContain(`height: 800px`)
  })

  it.each(Object.values(MD_METRICS))(
    `$key label points at metrics.md and has correct direction`,
    (label) => {
      expect(label.path).toBe(`metrics.md`)
      // combined_score (CMDS) is a score (higher=better); the rest are errors (lower)
      const expected = label.key === `combined_score` ? `higher` : `lower`
      expect(label.better).toBe(expected)
    },
  )

  it(`restores URL state for scatter axes and table sort`, async () => {
    await mount_with_url(
      MdPage,
      `http://localhost/tasks/md?x=combined_score&y=force_rmse&sort=rdf_error&dir=asc`,
    )

    const heading = document.querySelector(`h2`)?.textContent?.replaceAll(/\s+/g, ` `)
    expect(heading).toContain(`FRMSE vs CMDS`)
    const header = sorted_header()
    expect(header?.textContent).toContain(`RDF`)
    expect(header?.getAttribute(`aria-sort`)).toBe(`ascending`)
  })
})
