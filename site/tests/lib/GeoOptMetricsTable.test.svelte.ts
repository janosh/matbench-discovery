import GeoOptMetricsTable from '$lib/GeoOptMetricsTable.svelte'
import { GEO_OPT_SYMMETRY_METRICS, HYPERPARAMS } from '$lib/labels'
import { mount, tick } from 'svelte'
import { describe, expect, it } from 'vitest'
import { doc_query } from '../index'

describe(`GeoOptMetricsTable`, () => {
  it(`renders table with correct structure, columns, groups, and units`, async () => {
    mount(GeoOptMetricsTable, { target: document.body })
    await tick()

    const table = document.querySelector(`table`)
    expect(table).not.toBeNull()
    expect(table?.querySelector(`thead`)).not.toBeNull()
    expect(table?.querySelector(`tbody`)).not.toBeNull()

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
    const symmetry_labels = Object.values(GEO_OPT_SYMMETRY_METRICS).map((m) => m.label)
    const found_symmetry = symmetry_labels.filter((label) =>
      header_html.some((html) => html.includes(label)),
    ).length
    expect(found_symmetry).toBeGreaterThan(0)

    // Group headers for Symmetry and Hyperparams
    const group_texts = [...document.querySelectorAll(`tr.group-header th`)]
      .map((h) => h.textContent?.trim())
      .filter(Boolean)
    expect(group_texts).toContain(`Symmetry`)
    expect(group_texts).toContain(`Hyperparams`)
  })

  it(`sets initial sort to RMSD ascending`, async () => {
    mount(GeoOptMetricsTable, {
      target: document.body,
      props: { show_compliant: true, show_non_compliant: true },
    })
    await tick()

    const headers = [...document.querySelectorAll(`th`)]
    const rmsd_header = headers.find((h) => h.innerHTML?.includes(`RMSD`))
    expect(rmsd_header).toBeDefined()
    expect(rmsd_header?.getAttribute(`aria-sort`)).toBe(`ascending`)
  })

  it(`renders rows with model data and links`, async () => {
    mount(GeoOptMetricsTable, {
      target: document.body,
      props: { show_compliant: true, show_non_compliant: true },
    })
    await tick()

    const rows = document.querySelectorAll(`tbody tr`)
    expect(rows.length).toBeGreaterThan(0)

    // Model cells should have links to model pages
    const model_cells = [...document.querySelectorAll(`td[data-col="Model"]`)]
    expect(model_cells.length).toBeGreaterThan(0)
    expect(model_cells[0]?.querySelector(`a`)?.getAttribute(`href`)).toMatch(
      /^\/models\//,
    )

    // RMSD cells should have numeric values (data-col includes enriched label with unit)
    const rmsd_cells = [...document.querySelectorAll(`td`)].filter((td) =>
      td.getAttribute(`data-col`)?.includes(`RMSD`),
    )
    expect(
      rmsd_cells.some((cell) => {
        const text = cell.textContent?.trim()
        return text && !isNaN(parseFloat(text))
      }),
    ).toBe(true)
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

  // TODO: re-enable after matterviz release exports ToggleMenu with testable DOM
  it.skip(`opens column visibility panel with checkboxes`, async () => {
    mount(GeoOptMetricsTable, { target: document.body })
    await tick()

    const toggle_btn = doc_query<HTMLElement>(`.column-toggles summary`)
    toggle_btn.click()
    await tick()

    const column_menu = document.querySelector(`.column-menu`)
    expect(column_menu).not.toBeNull()
    expect(
      column_menu?.querySelectorAll(`input[type="checkbox"]`).length,
    ).toBeGreaterThan(0)
  })

  it.each([`RMSD`, `Model`])(`sorts by %s when header is clicked`, async (col_name) => {
    mount(GeoOptMetricsTable, {
      target: document.body,
      props: { show_compliant: true, show_non_compliant: true },
    })
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

    if (initial.length > 1) {
      const changed =
        JSON.stringify(initial) !== JSON.stringify(after_click) ||
        JSON.stringify(after_click) !== JSON.stringify(after_second)
      expect(changed).toBe(true)
    }
  })

  it.each([
    { show_compliant: true, show_non_compliant: true, desc: `all models` },
    { show_compliant: true, show_non_compliant: false, desc: `compliant only` },
    { show_compliant: false, show_non_compliant: true, desc: `non-compliant only` },
    { show_compliant: false, show_non_compliant: false, desc: `no models` },
  ])(`renders correctly with $desc`, async ({ show_compliant, show_non_compliant }) => {
    mount(GeoOptMetricsTable, {
      target: document.body,
      props: { show_compliant, show_non_compliant },
    })
    await tick()

    const table = document.querySelector(`table`)
    expect(table).not.toBeNull()

    const rows = document.querySelectorAll(`tbody tr`)
    if (!show_compliant && !show_non_compliant) {
      // HeatmapTable may render a "no data" placeholder row when empty
      expect(rows.length).toBeLessThanOrEqual(1)
    }
  })

  it.each([
    { prop: `show_compliant`, label_match: /^Compliant$|Compliant(?! only)/ },
    { prop: `show_non_compliant`, label_match: /Non-compliant/ },
  ])(`binds $prop prop correctly`, async ({ prop, label_match }) => {
    const state: Record<string, boolean> = { [prop]: true }
    mount(GeoOptMetricsTable, {
      target: document.body,
      props: {
        get [prop]() {
          return state[prop]
        },
        set [prop](val: boolean) {
          state[prop] = val
        },
      },
    })
    await tick()

    const labels = [...document.querySelectorAll(`label`)]
    const label = labels.find((l) => label_match.test(l.textContent ?? ``))
    const checkbox = label?.querySelector<HTMLInputElement>(`input[type="checkbox"]`)

    expect(state[prop]).toBe(true)
    if (checkbox) {
      checkbox.click()
      await tick()
      expect(state[prop]).toBe(false)
    }
  })
})
