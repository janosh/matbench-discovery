import { ACTIVE_MODELS } from '$lib'
import { MD_METRICS } from '$lib/labels'
import MdPage from '$routes/tasks/md/+page.svelte'
import { describe, expect, it } from 'vitest'
import {
  checkbox_for,
  doc_query,
  filter_summary_badge,
  has_md_metrics,
  mount,
  mount_with_url,
  sorted_header,
} from '../index'

describe(`MD Task Page`, () => {
  it(`renders page structure with filtered leaderboard and scatter`, () => {
    mount(MdPage, { target: document.body })

    expect(document.querySelector(`h1`)?.textContent).toContain(
      `Molecular Dynamics Metrics`,
    )

    const table = doc_query(`section.full-bleed table`)
    expect(table.querySelectorAll(`tbody tr`)).toHaveLength(
      ACTIVE_MODELS.filter(has_md_metrics).length,
    )
    const headers = [...table.querySelectorAll(`th`)].map((th) =>
      th.textContent?.replace(/\s*[↑↓]\s*$/, ``).trim(),
    )
    const expected = [`ΔERMSE`, `FRMSE`, `ΔADF`, `ΔvDOS`, `PMAE`, `PW1`, `CMDS`]
    for (const header of [...expected, `Speed`, `Slowdown`]) {
      expect(headers, `missing column ${header}`).toContain(header)
    }
    // ΔRDF is hidden from the leaderboard (redundant with ΔvDOS/ΔADF, out of CMDS);
    // the diatomics Speed/Slowdown columns must not leak onto the MD page despite
    // sharing labels with the MD ones (col visibility is keyed by unique col.key)
    expect(headers).not.toContain(`ΔRDF`)
    expect(headers.filter((header) => header === `Speed`)).toHaveLength(1)
    expect(headers.filter((header) => header === `Slowdown`)).toHaveLength(1)

    const headings = [...document.querySelectorAll<HTMLHeadingElement>(`h2`)].map((h2) =>
      h2.textContent?.replaceAll(/\s+/g, ` `).trim(),
    )
    expect(headings).toContain(`CMDS vs Speed`)

    const scatter = doc_query<HTMLDivElement>(`div.scatter`)
    expect(scatter.getAttribute(`style`)).toContain(`height: 800px`)
  })

  it.each(Object.values(MD_METRICS))(
    `$key label points at metrics.md and has correct direction`,
    (label) => {
      // Slowdown is computed per table view in assemble_row_data, not read from models
      const expected_path = label.key === `md_time_multiplier` ? undefined : `metrics.md`
      expect(label.path).toBe(expected_path)
      // combined_score (CMDS) is a score (higher=better); the rest are errors (lower)
      const expected = label.key === `combined_score` ? `higher` : `lower`
      expect(label.better).toBe(expected)
    },
  )

  it(`restores URL state for scatter axes, table sort and filters`, async () => {
    await mount_with_url(
      MdPage,
      `http://localhost/tasks/md?x=combined_score&y=force_rmse&sort=vdos_error&dir=asc&train=-OMat24&heatmap=0`,
    )

    const heading = document.querySelector(`h2`)?.textContent?.replaceAll(/\s+/g, ` `)
    expect(heading).toContain(`FRMSE vs CMDS`)
    const header = sorted_header()
    expect(header?.textContent).toContain(`vDOS`)
    expect(header?.getAttribute(`aria-sort`)).toBe(`ascending`)
    expect(filter_summary_badge(`Training data`)).toContain(`(1)`)
    expect(checkbox_for(`Heatmap`).checked).toBe(false)
  })
})
