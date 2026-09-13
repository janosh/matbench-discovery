import { goto } from '$app/navigation'
import DATASETS from '$data/datasets.yml'
import { ACTIVE_MODELS, make_table_filters } from '$lib/models.svelte'
import Page from '$routes/data/sets/+page.svelte'
import type { ScatterPlot } from 'matterviz/plot'
import { tick, type ComponentProps } from 'svelte'
import { describe, expect, it, vi } from 'vitest'
import {
  choose_scatter_property,
  doc_query,
  mount,
  mount_with_url,
  sorted_header,
} from '../index'

const plot_mock = vi.hoisted(() => vi.fn())
vi.mock(`matterviz/plot`, async (import_original) => ({
  ...(await import_original<Record<string, unknown>>()),
  ScatterPlot: plot_mock,
}))

const rows = () => [
  ...document.querySelectorAll<HTMLTableRowElement>(`.heatmap tbody tr[data-row-idx]`),
]
const names = () =>
  rows().map((row) => row.querySelector(`[data-col="Name"] > a`)?.textContent)
const dataset_row = (key: string) => {
  const row = doc_query(
    `.heatmap [data-col="Name"] > a[href="/data/${DATASETS[key].slug}"]`,
  ).closest(`tr`)
  if (!row) throw new Error(`Missing dataset row: ${key}`)
  return row
}
const dataset_cell = (key: string, column: string) =>
  doc_query(`[data-col="${column}"]`, dataset_row(key))
const plot_props = () =>
  plot_mock.mock.calls.at(-1)?.[1] as ComponentProps<typeof ScatterPlot>

describe(`Datasets Page`, () => {
  it(`renders all datasets with informative metadata and a permanent plot`, () => {
    plot_mock.mockClear()
    mount(Page, { target: document.body })
    expect(rows()).toHaveLength(Object.keys(DATASETS).length)
    expect(document.querySelectorAll(`.heatmap th`)).toHaveLength(12)
    expect(sorted_header()?.textContent).toContain(`Models`)
    expect(sorted_header()?.getAttribute(`aria-sort`)).toBe(`descending`)
    expect(names()[0]).toBe(`MPtrj`)
    const model_counts = rows().map((row) =>
      Number(doc_query(`[data-col="Models"]`, row).textContent),
    )
    expect(model_counts).toEqual(model_counts.toSorted((left, right) => right - left))
    expect(doc_query(`#dataset-growth h2`).textContent).toBe(`Dataset Sizes Over Time`)
    expect(document.querySelector(`details#dataset-growth`)).toBeNull()
    expect(
      document.querySelectorAll(`.property-picker input[role="combobox"]`),
    ).toHaveLength(1)
    for (const axis of [plot_props().x_axis, plot_props().y_axis]) {
      expect(axis?.options?.map(({ label }) => label)).toEqual([
        `Created`,
        `Structures`,
        `Materials`,
        `Models`,
      ])
    }
    expect(plot_props().series?.flatMap(({ x }) => x)).toHaveLength(
      Object.values(DATASETS).filter(({ n_structures }) => n_structures !== null).length,
    )
    expect(
      plot_props()
        .series?.flatMap(({ y }) => y)
        .every(Number.isFinite),
    ).toBe(true)
    for (const [key, role, access, release] of [
      [`WBM`, `Test`, `Public`, `Fixed release`],
      [`sAlex Validation`, `Validation`, `Public`, `Fixed release`],
      [`NOMAD`, `Repository`, `Public`, `Updated`],
      [`ELEMENTA`, `Training`, `Partial`, `Fixed release`],
      [`GNoME`, `Training`, `Partial`, `Fixed release`],
      [`MatterSim`, `Training`, `Unreleased`, `Fixed release`],
      [`MDR-MP PBE ω_q`, `Training`, `Public`, `Fixed release`],
      [`SMAX`, `Training`, `Unreleased`, `Fixed release`],
      [`OC25`, `Training`, `Public`, `Fixed release`],
      [`OMol25`, `Training`, `Public`, `Fixed release`],
      [`OMC25`, `Training`, `Public`, `Fixed release`],
      [`ODAC23`, `Training`, `Public`, `Fixed release`],
      [`ODAC25`, `Training`, `Public`, `Fixed release`],
    ]) {
      const row = dataset_row(key)
      for (const [column, value] of [
        [`Role`, role],
        [`Access`, access],
        [`Release`, release],
      ])
        expect(doc_query(`[data-col="${column}"]`, row).textContent?.trim()).toBe(value)
    }
    for (const [key, method, ending] of [
      [`MPtrj`, `DFT · PBE+U`, `PBE+U`],
      [`MatPES r2SCAN`, `DFT · r2SCAN`, `r2SCAN`],
      [`NOMAD`, `DFT, ML, GW, DMFT, MD · Various`, ` Various`],
      [`OC20`, `DFT · RPBE`, ` RPBE`],
      [`OC22`, `DFT · PBE+U`, `PBE+U`],
      [`OMol25`, `DFT · wB97M-V`, `B97M-V`],
      [`OpenLAM`, `DFT · Various`, `Various`.slice(1)],
    ]) {
      const cell = dataset_cell(key, `Method`)
      expect(cell.textContent?.trim()).toBe(method)
      expect(cell.style.maxWidth).toBe(`12em`)
      expect(doc_query(`.method`, cell).getAttribute(`title`)).toContain(method)
      expect(doc_query(`.method span:last-child`, cell).textContent).toBe(ending)
    }
    expect(dataset_cell(`MPtrj`, `Structures`).textContent).toMatch(/1\.58M/)
    expect(
      doc_query(`span`, dataset_cell(`SMAX`, `Materials`)).getAttribute(`title`),
    ).toContain(`not reported`)
    for (const [key, column, count] of [
      [`MAD-1.6`, `Structures`, `362,646`],
      [`OC20NEB`, `Structures`, `1,375,930`],
      [`OC20NEB`, `Materials`, `821`],
      [`MatPES PBE`, `Materials`, `110,317`],
      [`MatPES r2SCAN`, `Materials`, `89,719`],
      [`QMOF`, `Materials`, `20,372`],
      [`OMC25`, `Materials`, `218,841`],
      [`rMD17`, `Materials`, `10`],
      [`OCx24`, `Materials`, `19,406`],
    ]) {
      const cell = doc_query(`span`, dataset_cell(key, column))
      expect(cell.textContent).not.toBe(`n/a`)
      expect(cell.getAttribute(`title`)).toBe(count)
    }
    const unknown_count = doc_query(`span`, dataset_cell(`OCx24`, `Structures`))
    expect(unknown_count.textContent).toBe(`n/a`)
    expect(unknown_count.getAttribute(`title`)).toContain(`Structure count not reported`)
    expect(DATASETS.OCx24.notes?.Size).toMatch(/\S/)
    expect(
      [...document.querySelectorAll(`.heatmap th.not-sortable`)].map((header) =>
        header.textContent?.trim(),
      ),
    ).toEqual([`API`, `Links`])
  })

  it.each([
    [`q=phonon`, [`MDR-MP PBE ω_q`, `MDR Phonon PBE`]],
    [`access=partial`, [`GNoME`, `ELEMENTA`]],
    [`role=test`, [`WBM`, `MDR Phonon PBE`, `OCx24`, `rMD17`]],
    [`access=partial&sort=Created`, [`ELEMENTA`, `GNoME`]],
    [`role=repository&method=ML`, [`NOMAD`]],
    [
      `q=OMol25&access=public&method=DFT`,
      [`COSMOSDataset`, `OMol25`, `OMol25 Electronic`, `OPoly26`],
    ],
    [`q=%20PROJECT%20trajectories&access=public&role=training&method=DFT`, [`MPtrj`]],
  ])(`restores intersecting filters: %s`, async (params, expected) => {
    await mount_with_url(Page, `http://localhost/data/sets?${params}`)
    expect(names()).toEqual(expected)
    expect(doc_query(`[role="status"]`).textContent).toContain(
      `${expected.length} of ${Object.keys(DATASETS).length}`,
    )
  })

  it(`keeps unknown-size datasets searchable without plotting invented counts`, async () => {
    plot_mock.mockClear()
    await mount_with_url(Page, `http://localhost/data/sets?q=OCx24`)
    expect(names()).toEqual([`OCx24`])
    expect(plot_props().series).toEqual([])
    expect(document.body.textContent).toContain(
      `No datasets have values for all selected columns`,
    )
    plot_props().on_axis_change?.(`y`, `n_models`)
    await tick()
    expect(plot_props().series?.flatMap(({ y }) => y)).toEqual([0])
  })

  it(`syncs controls and sort, handles empty results, and resets only filters`, async () => {
    await mount_with_url(
      Page,
      `http://localhost/data/sets?access=invalid&role=invalid&method=invalid&sort=Links&dir=bad&color=invalid&x=role&y=invalid&size=access`,
    )
    expect(location.search).toBe(``)
    const search = doc_query<HTMLInputElement>(`[aria-label="Search datasets"]`)
    search.value = `no-such-dataset`
    search.dispatchEvent(new Event(`input`, { bubbles: true }))
    await tick()
    expect(rows()).toHaveLength(0)
    expect(document.body.textContent).toContain(`No datasets match`)
    expect(new URL(location.href).searchParams.get(`q`)).toBe(`no-such-dataset`)
    doc_query<HTMLButtonElement>(`.empty button`).click()
    await tick()
    expect(rows()).toHaveLength(Object.keys(DATASETS).length)
    for (const [label, key, value] of [
      [`Access`, `access`, `partial`],
      [`Role`, `role`, `training`],
      [`Method`, `method`, `DFT`],
    ]) {
      const select = doc_query<HTMLSelectElement>(`select[aria-label="${label}"]`)
      select.value = value
      // happy-dom does not match option:checked in Svelte's change handler.
      vi.spyOn(select, `querySelector`).mockReturnValueOnce(select.selectedOptions[0])
      select.dispatchEvent(new Event(`change`, { bubbles: true }))
      await tick()
      expect(new URLSearchParams(location.search).get(key)).toBe(value)
    }
    plot_props().color_bar?.on_property_change?.(`role`)
    await tick()
    expect(new URLSearchParams(location.search).get(`color`)).toBe(`role`)
    const header = [...document.querySelectorAll<HTMLTableCellElement>(`th`)].find(
      (candidate) => candidate.textContent?.trim() === `Name`,
    )
    header?.click()
    await tick()
    header?.click()
    await tick()
    expect(new URL(location.href).searchParams.get(`sort`)).toBe(`Name`)
    search.value = `MPtrj`
    search.dispatchEvent(new Event(`input`, { bubbles: true }))
    await tick()
    doc_query<HTMLButtonElement>(`.filters button`).click()
    await tick()
    expect(names()[0]).toBe(`AFLOW`)
    expect(new URL(location.href).searchParams.get(`sort`)).toBe(`Name`)
    for (const key of [`q`, `access`, `role`, `method`]) {
      expect(new URLSearchParams(location.search).has(key)).toBe(false)
    }
    expect(new URLSearchParams(location.search).get(`color`)).toBe(`role`)
    header?.click()
    await tick()
    expect(new URL(location.href).searchParams.get(`sort`)).toBe(`Name`)
    expect(sorted_header()?.getAttribute(`aria-sort`)).toBe(`descending`)
  })

  it(`links exact model counts to matching leaderboard filters without expanding composites`, () => {
    mount(Page, { target: document.body })
    for (const key of Object.keys(DATASETS)) {
      const count = ACTIVE_MODELS.filter((model) =>
        model.training_sets.includes(key),
      ).length
      const cell = dataset_cell(key, `Models`)
      expect(Number(cell.textContent)).toBe(count)
      const link = cell.querySelector(`a`)
      if (!count) {
        expect(link).toBeNull()
        continue
      }
      const filters = make_table_filters()
      filters.read(new URL(link?.href ?? ``).searchParams)
      expect(ACTIVE_MODELS.filter(filters.matches)).toHaveLength(count)
    }
    expect(Number(dataset_cell(`Alex`, `Models`).textContent)).toBeLessThan(
      ACTIVE_MODELS.filter((model) =>
        model.training_sets.some(
          (key) => key === `Alex` || DATASETS[key].contains?.includes(`Alex`),
        ),
      ).length,
    )
  })

  it(`opens explanations by click, focus, and hover with source and caveat links`, async () => {
    mount(Page, { target: document.body })
    for (const interaction of [`click`, `focus`, `hover`]) {
      const button = doc_query<HTMLButtonElement>(`button[aria-labelledby="about-salex"]`)
      if (interaction === `click`) button.click()
      else if (interaction === `focus`) button.focus()
      else button.dispatchEvent(new MouseEvent(`mouseenter`))
      await vi.waitFor(() => expect(button.getAttribute(`aria-expanded`)).toBe(`true`))
      // A direct aria-label would trigger a second tooltip from the table delegation.
      expect(button.hasAttribute(`aria-label`)).toBe(false)
      expect(doc_query(`#about-salex`).textContent).toBe(`About sAlex`)
      const content = doc_query(`[id="${button.getAttribute(`aria-controls`)}"]`)
      expect(content.textContent).toContain(`Derived from:`)
      expect(content.querySelector(`a[href="/data/alex"]`)).not.toBeNull()
      button.blur()
      button.dispatchEvent(new MouseEvent(`mouseleave`))
      await vi.waitFor(() => expect(button.getAttribute(`aria-expanded`)).toBe(`false`))
    }
    const button = doc_query<HTMLButtonElement>(`button[aria-labelledby="about-mptrj"]`)
    button.focus()
    await vi.waitFor(() => expect(button.getAttribute(`aria-expanded`)).toBe(`true`))
    expect(
      doc_query(`[id="${button.getAttribute(`aria-controls`)}"]`).textContent,
    ).toContain(`Leakage`)
  })

  it(`keeps resource and API links accessible and dataset-specific`, () => {
    mount(Page, { target: document.body })
    for (const [key, dataset] of Object.entries(DATASETS)) {
      const links = [
        ...dataset_row(key).querySelectorAll(`[data-col="API"] a, [data-col="Links"] a`),
      ]
      expect(
        links
          .map((link) => link.getAttribute(`href`))
          .toSorted((left, right) => String(left).localeCompare(String(right))),
      ).toEqual(
        [
          dataset.native_api,
          dataset.optimade_api,
          dataset.url,
          dataset.download_url,
          dataset.doi,
        ]
          .filter(Boolean)
          .toSorted((left, right) => String(left).localeCompare(String(right))),
      )
      for (const link of links) {
        expect(link.getAttribute(`target`)).toBe(`_blank`)
        expect(link.getAttribute(`rel`)).toContain(`noopener`)
        expect(link.getAttribute(`aria-label`)).toBe(
          `${link.getAttribute(`title`)} for ${key}`,
        )
      }
    }
  })

  it.each([
    [`access`, `Partial`, `#b16c00`],
    [`role`, `Training`, `#397ec5`],
  ])(
    `colors filtered datasets by %s and opens selected datasets`,
    async (color_by, label, color) => {
      await mount_with_url(
        Page,
        `http://localhost/data/sets?access=partial&color=${color_by}`,
      )
      const props = plot_props()
      expect(props).toMatchObject({
        x_axis: { scale_type: `time` },
        y_axis: { label: `Structures` },
        color_bar: { categories: { [label]: color }, selected_property_key: color_by },
        series: [`GNoME`, `ELEMENTA`].map((key) => ({
          label: key,
          x: [new Date(DATASETS[key].date_created).getTime()],
          y: [DATASETS[key].n_structures],
          point_style: { fill: color },
        })),
      })
      expect(props.series?.[0].point_label).toEqual({
        text: `GNoME`,
        font_size: `12px`,
        auto_placement: true,
      })
      const click = plot_props().point_events?.onclick
      const metadata = props.series?.[1].metadata
      if (!metadata || Array.isArray(metadata))
        throw new Error(`Expected dataset metadata`)
      expect(metadata).toMatchObject({ key: `ELEMENTA`, href: `/data/elementa` })
      click?.({
        point: { x: 0, y: 1, series_idx: 1, point_idx: 0, metadata },
        event: new MouseEvent(`click`),
      })
      expect(goto).toHaveBeenCalledWith(`/data/elementa`)
    },
  )

  it.each([
    [`x`, `X axis`, `x`],
    [`y`, `Y axis`, `y`],
    [`color`, `Color`, `color_values`],
    [`size`, `Marker size`, `size_values`],
  ] as const)(`selects and restores the %s column`, async (dim, label, series_key) => {
    await mount_with_url(
      Page,
      `http://localhost/data/sets?x=n_models&y=date_created&color=n_structures&size=n_models`,
    )
    expect(plot_props().x_axis?.label).toBe(`Models`)
    expect(plot_props().y_axis).toMatchObject({
      label: `Created`,
      scale_type: `time`,
      format: `%Y`,
    })
    expect(plot_props().color_bar).toMatchObject({ title: `Structures` })
    if (dim === `x` || dim === `y`) {
      plot_props().on_axis_change?.(dim, `n_materials`)
      await tick()
      expect(plot_props()[`${dim}_axis`]?.selected_key).toBe(`n_materials`)
    } else if (dim === `color`) {
      plot_props().color_bar?.on_property_change?.(`n_materials`)
      await tick()
      expect(plot_props().color_bar?.selected_property_key).toBe(`n_materials`)
    } else await choose_scatter_property(label, `Materials`)
    expect(new URLSearchParams(location.search).get(dim)).toBe(`n_materials`)
    const series = plot_props().series ?? []
    const entries = Object.entries(DATASETS).filter(
      ([, dataset]) =>
        dataset.n_materials != null && (dim === `color` || dataset.n_structures !== null),
    )
    expect(series.map(({ id }) => id)).toEqual(entries.map(([key]) => key))
    expect(series.flatMap((entry) => entry[series_key])).toEqual(
      entries.map(([, dataset]) => dataset.n_materials),
    )
    if (dim === `color`) expect(plot_props().color_bar?.categories).toBeUndefined()
  })
})
