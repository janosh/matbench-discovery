import { ACTIVE_MODELS } from '$lib'
import GeoOptMetricsTable from '$lib/table/GeoOptMetricsTable.svelte'
import { GEO_OPT_SYMMETRY_METRICS, HYPERPARAMS } from '$lib/labels'
import { make_table_filters } from '$lib/models.svelte'
import type { ModelData } from '$lib/types'
import { tick } from 'svelte'
import { describe, expect, it } from 'vitest'
import { doc_query, mount } from '../index'

// Mirrors the component's geo-opt presence check and default table filters.
const geo_opt_row_count = (
  matches: (model: ModelData) => boolean = make_table_filters().matches,
) =>
  ACTIVE_MODELS.filter((model) => model.metrics?.geo_opt != null && matches(model)).length

describe(`GeoOptMetricsTable`, () => {
  it(`renders table with correct structure, columns, groups, and units`, async () => {
    mount(GeoOptMetricsTable, { target: document.body })
    await tick()

    const table = doc_query(`table`)
    doc_query(`thead`, table)
    doc_query(`tbody`, table)

    const headers = [...document.querySelectorAll(`th`)]
    const header_texts = headers.map((h) => h.textContent?.trim())
    const header_html = headers.map((h) => h.innerHTML)

    // Model column present
    expect(header_texts).toContain(`Model`)

    // RMSD header omits the unit for a concise default column header
    const rmsd_header = header_html.find((h) => h.includes(`RMSD`))
    expect(rmsd_header).not.toContain(`unitless`)
    expect(rmsd_header).not.toContain(`font-weight: 200`)

    // f_max header includes unit
    const fmax_header = header_html.find((h) => h.includes(`f<sub>max</sub>`))
    expect(fmax_header).toContain(`(${HYPERPARAMS.max_force.unit})`)

    // Visible hyperparams present
    expect(header_texts).toContain(HYPERPARAMS.ase_optimizer.label)
    expect(header_texts).toContain(HYPERPARAMS.max_steps.label)
    expect(header_texts).toContain(HYPERPARAMS.cell_filter.label)

    // Hidden-by-default columns absent
    expect(header_texts).not.toContain(HYPERPARAMS.n_layers.label)
    expect(
      header_html.some((h) => h.includes(HYPERPARAMS.graph_construction_radius.label)),
    ).toBe(false)

    // At least some symmetry metrics present
    const symmetry_labels = Object.values(GEO_OPT_SYMMETRY_METRICS).map(
      (metric) => metric.label,
    )
    const found_symmetry = symmetry_labels.filter((label) =>
      header_html.some((html) => html.includes(label)),
    ).length
    expect(found_symmetry).toBe(symmetry_labels.length)

    // Group headers for Symmetry and Hyperparams
    const group_texts = [...document.querySelectorAll(`tr.group-header th`)]
      .map((h) => h.textContent?.trim())
      .filter(Boolean)
    expect(group_texts).toContain(`Symmetry`)
    expect(group_texts).toContain(`Hyperparams`)
  })

  it(`sets initial sort to RMSD ascending`, async () => {
    mount(GeoOptMetricsTable, { target: document.body })
    await tick()

    const rmsd_header = [...document.querySelectorAll(`th`)].find((header) =>
      header.innerHTML.includes(`RMSD`),
    )
    if (!rmsd_header) throw new Error(`RMSD header not found`)
    expect(rmsd_header.getAttribute(`aria-sort`)).toBe(`ascending`)
  })

  it(`renders geo-opt rows and excludes models without geo-opt metrics`, async () => {
    const model_key = `no-geo-opt-regression`
    ACTIVE_MODELS.push({
      model_key,
      model_name: model_key,
      model_version: `test`,
      targets: `EFS_G`,
      training_sets: [],
      n_training_materials: 1,
      n_training_structures: 1,
      model_params: 1,
      dates: { benchmark_added: `2026-06-30` },
      metrics: {
        discovery: { full_test_set: { F1: 0.1 } },
        diatomics: { energy_mae: 1 },
      },
    } as unknown as ModelData)

    try {
      mount(GeoOptMetricsTable, { target: document.body })
      await tick()

      expect(document.body.textContent).not.toContain(model_key)

      const rows = document.querySelectorAll(`tbody tr`)
      expect(rows).toHaveLength(geo_opt_row_count())

      const model_cells = [...document.querySelectorAll(`td[data-col="Model"]`)]
      expect(model_cells).toHaveLength(rows.length)
      expect(model_cells[0]?.querySelector(`a`)?.getAttribute(`href`)).toMatch(
        /^\/models\//,
      )

      const rmsd_cells = [...document.querySelectorAll(`td`)].filter((td) =>
        td.getAttribute(`data-col`)?.includes(`RMSD`),
      )
      expect(
        rmsd_cells.some((cell) => {
          const text = cell.textContent?.trim()
          return text ? !Number.isNaN(Number(text)) : false
        }),
      ).toBe(true)
    } finally {
      ACTIVE_MODELS.pop()
    }
  })

  it(`toggles heatmap colors`, async () => {
    mount(GeoOptMetricsTable, { target: document.body })
    await tick()

    const checkbox = doc_query<HTMLInputElement>(
      `input[type="checkbox"][aria-label="Toggle heatmap colors"]`,
    )
    expect(checkbox.checked).toBe(true)

    checkbox.click()
    await tick()
    expect(checkbox.checked).toBe(false)

    checkbox.click()
    await tick()
    expect(checkbox.checked).toBe(true)
  })

  it(`initializes and binds column_order`, async () => {
    const state = { column_order: [] as string[] }
    mount(GeoOptMetricsTable, {
      target: document.body,
      props: {
        get column_order() {
          return state.column_order
        },
        set column_order(val) {
          state.column_order = val
        },
      },
    })
    await tick()

    expect(state.column_order.length).toBeGreaterThan(5)
    expect(state.column_order).toContain(`Model`)
  })

  it(`opens column visibility panel with checkboxes`, async () => {
    mount(GeoOptMetricsTable, { target: document.body })
    await tick()

    const toggle_btn = doc_query(`.column-toggles summary`)
    toggle_btn.click()
    await tick()

    const column_menu = doc_query(`.sections-container`)
    expect(
      [...column_menu.querySelectorAll(`.section-header`)].map((button) =>
        button.textContent?.trim().replace(/^[▶▼]\s*/, ``),
      ),
    ).toStrictEqual([`Symmetry`, `Hyperparams`])
    expect(column_menu.querySelectorAll(`input[type="checkbox"]`).length).toBeGreaterThan(
      Object.keys(GEO_OPT_SYMMETRY_METRICS).length,
    )
  })

  it.each([`RMSD`, `Model`])(`sorts by %s when header is clicked`, async (col_name) => {
    mount(GeoOptMetricsTable, { target: document.body })
    await tick()

    const headers = [...document.querySelectorAll(`th`)]
    const header = headers.find((h) =>
      col_name === `Model`
        ? h.textContent?.trim() === `Model`
        : h.textContent?.includes(col_name),
    )
    if (!header) throw new Error(`${col_name} header not found`)

    const get_order = () =>
      [...document.querySelectorAll(`td[data-col="Model"]`)].map((c) =>
        c.textContent?.trim(),
      )

    const initial = get_order()
    header.click()
    await tick()
    const after_click = get_order()
    header.click()
    await tick()
    const after_second = get_order()

    const changed =
      JSON.stringify(initial) !== JSON.stringify(after_click) ||
      JSON.stringify(after_click) !== JSON.stringify(after_second)
    expect(initial).toHaveLength(geo_opt_row_count())
    expect(changed).toBe(true)
  })

  it.each([
    { training: {}, desc: `no filters` },
    { training: { MPtrj: `require` }, desc: `MPtrj-trained only` },
    { training: { OMat24: `exclude` }, desc: `OMat24 excluded` },
  ] as const)(`filters rows with $desc`, async ({ training }) => {
    const filters = make_table_filters()
    filters.training = { ...training } as typeof filters.training
    mount(GeoOptMetricsTable, { target: document.body, props: { filters } })
    await tick()

    doc_query(`thead`, doc_query(`table`))

    const rows = document.querySelectorAll(`tbody tr`)
    const expected_rows = geo_opt_row_count(filters.matches)
    // HeatmapTable may render a "no data" placeholder row when empty
    if (expected_rows === 0) expect(rows.length).toBeLessThanOrEqual(1)
    else expect(rows).toHaveLength(expected_rows)
  })
})
