import { ACTIVE_MODELS } from '#lib/models.svelte.js'
import { DEFAULT_CMDS_CONFIG } from '#lib/combined-scores.svelte.js'
import { is_finite_num } from '#lib/metrics.js'
import MdPage from '#routes/benchmarks/md/+page.svelte'
import { describe, expect, it } from 'vite-plus/test'
import {
  doc_query,
  expect_table_url_state,
  expect_task_page_layout,
  filter_summary_badge,
  has_md_metrics,
  header_name,
  heading_texts,
  mount,
  mount_with_url,
  patch_model,
  REPRO_CASES,
} from '../index.js'

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
    const headers = [...table.querySelectorAll(`th`)].map(header_name)
    const expected = [`ΔERMSE`, `FRMSE`, `ΔADF`, `ΔvDOS`, `PMAE`, `PW1`, `CMDS`]
    for (const header of [...expected, `twall`, `t/tmin`]) {
      expect(headers, `missing column ${header}`).toContain(header)
    }
    // ΔRDF is hidden from the leaderboard (redundant with ΔvDOS/ΔADF, out of CMDS);
    // the diatomics wall-time columns must not leak onto the MD page despite
    // sharing labels with the MD ones (col visibility is keyed by unique col.key)
    expect(headers).not.toContain(`ΔRDF`)
    expect(headers.filter((header) => header === `twall`)).toHaveLength(1)
    expect(headers.filter((header) => header === `t/tmin`)).toHaveLength(1)

    expect_task_page_layout(
      [`Leaderboard`, `CMDS vs twall`, `Test set: DynaMat v1.0`, `Methodology`],
      Object.keys(DEFAULT_CMDS_CONFIG).length,
    )
    const methodology = doc_query(`#methodology + p`)
    expect(table.compareDocumentPosition(methodology)).toBe(
      Node.DOCUMENT_POSITION_FOLLOWING,
    )
    expect(methodology.textContent).toContain(`same initial structures`)
    expect(doc_query(`a[href="#test-set"]`).compareDocumentPosition(table)).toBe(
      Node.DOCUMENT_POSITION_FOLLOWING,
    )

    const scatter = doc_query<HTMLDivElement>(`div.scatter`)
    expect(scatter.getAttribute(`style`)).toContain(`height: 800px`)
  })

  it(`restores URL state for scatter axes, table sort and filters`, async () => {
    await mount_with_url(
      MdPage,
      `http://localhost/benchmarks/md?x=combined_score&y=force_rmse&sort=vdos_error&dir=asc&train=-OMat24&heatmap=0`,
    )

    expect(heading_texts(`#model-comparison`)).toEqual([`FRMSE vs CMDS`])
    await expect_table_url_state(`ΔvDOS`, `ascending`)
    expect(filter_summary_badge(`Training data`)).toContain(`(1)`)
    // the scatter applies the table filters too (its size picker counts plotted models)
    const sized = ACTIVE_MODELS.filter(
      (model) => has_md_metrics(model) && is_finite_num(model.model_params),
    )
    const shown = sized.filter((model) => !model.training_sets.includes(`OMat24`))
    expect(shown.length).toBeLessThan(sized.length)
    expect(doc_query(`.property-picker .selected-label small`).textContent).toMatch(
      new RegExp(`\\b${shown.length}\\s+models`),
    )
  })

  it(`shows a Reproducible status column last`, () => {
    const model = ACTIVE_MODELS.find(has_md_metrics)
    if (!model) throw new Error(`No model with MD metrics`)
    const md = { ...model.metrics?.md, reproducibility: REPRO_CASES.drifted }
    patch_model(model, { metrics: { ...model.metrics, md } })
    mount(MdPage, { target: document.body })

    const headers = [...document.querySelectorAll(`section.full-bleed thead th`)]
    expect(header_name(headers.at(-1) ?? document.body)).toBe(`Reproducible`)
    const row = doc_query(`section.full-bleed a[href="/models/${model.model_key}"]`)
    expect(row.closest(`tr`)?.querySelector(`td:last-child`)?.textContent?.trim()).toBe(
      `≠ Drifted`,
    )
  })
})
