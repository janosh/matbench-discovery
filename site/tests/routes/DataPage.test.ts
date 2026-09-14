import Page from '$routes/data/data-files-direct-download.md'
import DataRoute from '$routes/data/+page.svelte'
import WbmDetails from '$routes/data/[slug]/WbmDetails.svelte'
import data_files from '$pkg/data-files.yml'
import TestSet from '$lib/benchmark/TestSet.svelte'
import { benchmarks } from '$lib/benchmark/data'
import { goto } from '$app/navigation'
import { tick } from 'svelte'
import { describe, expect, it, vi } from 'vitest'
import { doc_query, mount, mount_with_url } from '../index'

it.each([
  [``, `/benchmarks`],
  [`#wbm`, `/benchmarks/discovery#test-set`],
  [`#wbm-title`, `/benchmarks/discovery#test-set`],
  [`#phonondb`, `/benchmarks/phonons#test-set`],
  [`#dynamat`, `/benchmarks/md#test-set`],
  [`#diatomics`, `/benchmarks/diatomics#test-set`],
  [`#training-data`, `/data/sets`],
  [`#downloading-data`, `/benchmarks#downloading-data`],
  [`#constructor`, `/benchmarks#constructor`],
  [`#toString`, `/benchmarks#toString`],
  [`#__proto__`, `/benchmarks#__proto__`],
])(`redirects the former data page %s to %s`, async (hash, destination) => {
  await mount_with_url(DataRoute, `http://localhost/data?cps_weights=1,0,0${hash}`)
  const target = new URL(destination, `http://localhost`)
  target.search = `?cps_weights=1,0,0`
  expect(vi.mocked(goto)).toHaveBeenCalledWith(
    `${target.pathname}${target.search}${target.hash}`,
    { replaceState: true },
  )
})

it.each([
  [`discovery`, `O`, 27946, 85, `WBM reference hull-distance distribution`],
  [`geo-opt`, `O`, 27946, 85, `WBM structure sizes`],
  [`phonons`, `Te`, 11, 34, `PhononDB reference conductivity distribution`],
  [`md`, `H`, 6, 22, `DynaMat reference systems`],
  [`diatomics`, `H`, 1, 92, null],
] as const)(
  `shows %s test-set coverage, credits, downloads, and reference EDA`,
  async (task, symbol, count, n_elements, eda) => {
    vi.spyOn(HTMLElement.prototype, `clientWidth`, `get`).mockReturnValue(600)
    vi.spyOn(HTMLElement.prototype, `clientHeight`, `get`).mockReturnValue(300)
    mount(TestSet, { target: document.body, props: { task } })
    await tick()
    const { dataset } = benchmarks[task]
    const section = doc_query(`section[aria-labelledby="test-set"]`)
    expect(doc_query(`h2`, section).textContent).toBe(`Test set: ${dataset.name}`)
    expect(doc_query(`.coverage`, section).textContent).toBe(dataset.coverage)
    expect(doc_query(`.availability`, section).textContent).toBe(dataset.availability)
    if (dataset.credit) {
      const credit = doc_query(`.credit a`, section)
      expect(credit.textContent).toBe(dataset.credit[0])
      expect(credit.getAttribute(`href`)).toBe(dataset.credit[1])
    }
    expect(
      [...section.querySelectorAll(`.links a`)].map((link) => link.getAttribute(`href`)),
    ).toEqual(dataset.links.map(([, href]) => href))
    if (eda) expect(section.querySelector(`[aria-label="${eda}"]`)).not.toBeNull()
    if (task === `geo-opt`) {
      const plots = section.querySelectorAll(`.distributions .bar-plot`)
      expect(plots).toHaveLength(2)
      const ticks = [...plots[0].querySelectorAll(`.y-axis .tick text`)].map(
        (label) => label.textContent,
      )
      expect(ticks).toEqual(expect.arrayContaining([`1`, `10`, `100`, `1k`, `10k`]))
      expect(ticks).not.toContain(`0`)
      expect(doc_query(`.y-axis .tick text`, plots[1]).textContent).toBe(`0`)
    }
    if (task === `md`) {
      expect(section.querySelectorAll(`tbody tr`)).toHaveLength(17)
      expect(section.textContent).toContain(`omit energies and forces`)
    }
    if (task === `phonons`) {
      const bars = [...section.querySelectorAll(`.histogram-series path`)]
      expect(bars.length).toBeGreaterThan(15)
      const total = bars.reduce(
        (sum, bar) => sum + Number(bar.getAttribute(`aria-label`)?.split(`, count `)[1]),
        0,
      )
      expect(total).toBe(103)
    }
    const table = doc_query(`.periodic-table`, section)
    expect(table.querySelector(`.element-tile:not([data-element-symbol])`)).toBeNull()
    const tiles = [...table.querySelectorAll<HTMLElement>(`[data-element-symbol]`)]
    expect(tiles.filter((tile) => tile.style.opacity !== `0.15`)).toHaveLength(n_elements)
    const ticks = [...table.querySelectorAll(`.tick-label`)].map((label) =>
      Number(label.textContent),
    )
    expect(ticks[0]).toBe(1)
    expect(ticks.at(-1)).toBe(count)
    expect(ticks.every(Number.isInteger)).toBe(true)
    expect(new Set(ticks).size).toBe(ticks.length)
    if (task === `diatomics`)
      expect(doc_query(`.colorbar`, table).textContent).not.toContain(`log`)
    // Count occurrences per structure/system/pair, not atoms or trajectory frames.
    for (const [element, expected_count] of [
      [symbol, count],
      [`Og`, 0],
    ] as const) {
      const tile = doc_query(`[data-element-symbol="${element}"]`, table)
      tile.dispatchEvent(new MouseEvent(`mouseenter`))
      await tick()
      expect(doc_query(`.tooltip`, table).textContent).toContain(
        `${dataset.count_unit} containing ${element}: ${expected_count.toLocaleString(`en-US`)}`,
      )
      tile.dispatchEvent(new MouseEvent(`mouseleave`))
      await tick()
    }
  },
)

describe(`Public data downloads`, () => {
  it(`renders every public registry file with links and descriptions`, () => {
    mount(Page, { target: document.body })
    const data_files_list = doc_query(`ol.data-files-list`)
    const list_items = [...data_files_list.children]
    const list_text = data_files_list.textContent
    const public_files = Object.entries(data_files).filter(
      ([key]) => !key.startsWith(`_`),
    )

    expect(list_items).toHaveLength(public_files.length)
    for (const [idx, [key, file]] of public_files.entries()) {
      if (typeof file === `string`) throw new Error(`Expected a data file: ${key}`)
      expect(doc_query(`strong code`, list_items[idx]).textContent).toBe(key)
      expect(doc_query(`a`, list_items[idx]).getAttribute(`href`)).toBe(file.url)
      expect(doc_query(`p`, list_items[idx]).textContent).toMatch(/\S/)
    }
    expect(list_text).not.toContain(`_private_reference_data`)
    expect(list_text).not.toContain(`private_labeled_reference`)

    for (const link of data_files_list.querySelectorAll(`a`)) {
      const url = link.getAttribute(`href`)
      expect(url).toMatch(
        /^https:\/\/(?:figshare\.com|.*materialsproject\.(?:org|com)|github\.com)\/.*$/,
      )

      expect(link.textContent).toMatch(/\S/)
    }
  })
})

describe(`WBM details URL state`, () => {
  const count_mode_text = () =>
    doc_query(`.multiselect:has(#count-mode) ul.selected`).textContent

  it.each([
    [``, `occurrence`],
    [`?count_mode=composition&color_scale=interpolatePlasma`, `composition`],
    [`?count_mode=bogus`, `occurrence`], // invalid value falls back to default
  ])(`restores count mode from URL %s`, async (query, expected_mode) => {
    await mount_with_url(WbmDetails, `http://localhost/data/wbm${query}`)

    expect(document.querySelectorAll(`.sunburst`)).toHaveLength(2)
    expect(count_mode_text()).toContain(expected_mode)
    expect(document.querySelector(`.periodic-table .colorbar`)).not.toBeNull()
    expect(new URL(location.href).searchParams.get(`color_scale`)).toBe(
      new URLSearchParams(query).get(`color_scale`),
    )
    const scale_picker = doc_query(`.multiselect:has(input[aria-label="Color scale"])`)
    const expected_scale = (
      new URLSearchParams(query).get(`color_scale`) ?? `interpolateViridis`
    ).replace(`interpolate`, ``)
    expect(doc_query(`ul.selected`, scale_picker).textContent).toContain(expected_scale)
    doc_query(`input[aria-label="Color scale"]`).focus()
    await tick()
    const cividis_option = [
      ...scale_picker.querySelectorAll<HTMLElement>(`ul.options li[aria-posinset]`),
    ].find((option) => option.textContent?.includes(`Cividis`))
    if (!cividis_option) throw new Error(`Missing Cividis option`)
    expect(cividis_option.querySelector(`.colorbar`)).not.toBeNull()
    cividis_option.click()
    await tick()
    expect(new URL(location.href).searchParams.get(`color_scale`)).toBe(
      `interpolateCividis`,
    )

    const picker = doc_query(`.multiselect:has(#count-mode)`)
    expect(picker.querySelector(`button[title^="Remove"]`)).toBeNull()
    doc_query(`#count-mode`).click()
    await tick()
    const other_mode = expected_mode === `occurrence` ? `composition` : `occurrence`
    const option = [...picker.querySelectorAll<HTMLElement>(`li[role="option"]`)].find(
      (element) => element.textContent?.includes(other_mode),
    )
    if (!option) throw new Error(`Missing count mode ${other_mode}`)
    option.click()
    await tick()
    expect(count_mode_text()).toContain(other_mode)
    expect(count_mode_text()).not.toContain(expected_mode)
    expect(picker.querySelectorAll(`ul.selected > li`)).toHaveLength(1)
  })

  const log_checkboxes = (): HTMLInputElement[] => [
    ...document.querySelectorAll<HTMLInputElement>(`.table-inset input[type="checkbox"]`),
  ]

  // the three element-count heatmaps share one log toggle (two-way bound), also in the URL
  it(`shares the log toggle across all heatmaps and syncs it to the URL`, async () => {
    await mount_with_url(WbmDetails, `http://localhost/data/wbm`)
    const checkboxes = log_checkboxes()
    expect(checkboxes).toHaveLength(3)
    expect(checkboxes.map((box) => box.checked)).toEqual([false, false, false])

    checkboxes[1].click()
    await tick()
    expect(checkboxes.map((box) => box.checked)).toEqual([true, true, true])
    expect(new URL(location.href).searchParams.get(`log`)).toBe(`1`)

    document.body.innerHTML = ``
    await mount_with_url(WbmDetails, `http://localhost/data/wbm?log=1`)
    expect(log_checkboxes().map((box) => box.checked)).toEqual([true, true, true])
  })
})
