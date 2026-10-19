import kappa_103_analysis from '#figs/kappa-103-analysis.jsonl'
import { ACTIVE_MODELS, make_table_filters } from '#lib/models.svelte.js'
import type * as KappaParity from '#lib/parity/kappa-parity.js'
import PhononsPage from '#routes/benchmarks/phonons/+page.svelte'
import PhononFreqParity from '#routes/benchmarks/phonons/PhononFreqParity.svelte'
import { tick } from 'svelte'
import { afterEach, describe, expect, it, vi } from 'vite-plus/test'
import {
  doc_query,
  expect_table_url_state,
  filter_summary_badge,
  header_name,
  heading_texts,
  mock_element_size,
  mount,
  mount_with_url,
  sorted_header,
} from '../index.js'

// per-test override of the analysis payload's model list (null = real payload)
const analysis_mock = vi.hoisted(() => ({ models: null as unknown[] | null, error: `` }))
vi.mock(`#lib/parity/kappa-parity.js`, async (import_original) => {
  const actual = await import_original<typeof KappaParity>()
  return {
    ...actual,
    load_kappa_analysis: async () => {
      if (analysis_mock.error) throw new Error(analysis_mock.error)
      const analysis = await actual.load_kappa_analysis()
      return analysis_mock.models
        ? { ...analysis, models: analysis_mock.models }
        : analysis
    },
  }
})
afterEach(() => {
  analysis_mock.models = null
  analysis_mock.error = ``
})

const default_filters = make_table_filters()
// Mirrors the page's model_filter plus MetricsTable's default filters.
const phonon_leaderboard_count = ACTIVE_MODELS.filter(
  (model) => model.metrics?.phonons?.kappa_103 != null && default_filters.matches(model),
).length

const kappa_sort_select = (): HTMLSelectElement =>
  doc_query<HTMLSelectElement>(`.kappa-model-select select`)

// the model picker is a svelte-widgets MultiSelect; its selected chip lives in this ul
const kappa_selected_model = (): string | null =>
  doc_query(`.kappa-model-select ul[aria-label="selected options"]`).textContent

// Open the multiselect dropdown and return its options in listed order.
async function kappa_model_options(): Promise<HTMLElement[]> {
  doc_query<HTMLInputElement>(`.kappa-model-select input`).click()
  await tick()
  return Array.from(
    document.querySelectorAll<HTMLElement>(`.kappa-model-select ul[role="listbox"] > li`),
  )
}

describe(`Phonons Task Page`, () => {
  it(`renders the task narrative, scatter, and diagnostics`, async () => {
    mount(PhononsPage, { target: document.body })

    expect(doc_query(`h1`).textContent).toContain(`MLFF Phonon Modeling Metrics`)
    expect(document.body.textContent).toContain(`cancellation can hide errors`)
    expect(document.body.textContent).toMatch(/censored\s+to\s+2\s+because/)

    expect(heading_texts()).toEqual([
      `Leaderboard`,
      `Model Comparison: κSRME vs κSRE`,
      `Model Inspector`,
      `ML vs DFT Lattice Thermal Conductivity`,
      `Test set: PhononDB PBE`,
      `Methodology`,
    ])
    const methodology = doc_query(`#methodology + p`)
    expect(methodology.textContent).toContain(`PhononDB PBE test set`)
    expect(methodology.closest(`details`)).toBeNull()

    const scatter = doc_query<HTMLDivElement>(`div.scatter`)
    expect(scatter.getAttribute(`style`)).toContain(`height: 800px`)

    // SRME-vs-kappa scatter and frequency parity plot side by side.
    await vi.waitFor(
      () => expect(document.querySelector(`.diagnostics-grid`)).not.toBeNull(),
      { timeout: 10_000 },
    )
    expect(doc_query(`.diagnostics-grid`).querySelectorAll(`div.scatter`)).toHaveLength(2)
  })

  it(`keeps the frequency-parity metric annotation inside the plot area`, async () => {
    mock_element_size(800, 600)
    const [entry] = kappa_103_analysis.models
    mount(PhononFreqParity, {
      target: document.body,
      props: { entry, base: kappa_103_analysis },
    })
    const annotation = await vi.waitFor(() => doc_query(`.plot-annotation`))
    const annotation_frame = annotation.closest(`foreignObject`)
    const plot_area = doc_query(`clipPath rect`, annotation.closest(`svg`))
    // Keep the metric inside the data area when the legend reserves bottom padding.
    for (const attr of [`x`, `y`, `width`, `height`]) {
      expect(annotation_frame?.getAttribute(attr)).toBe(plot_area.getAttribute(attr))
    }
    expect(annotation.style.inset).toBe(`auto 0.5em 0.5em auto`)
  })

  // a leaderboard model absent from kappa-103-analysis.jsonl used to leave the spinner up forever
  it.each([false, true])(
    `distinguishes missing diagnostics from module failure (%s)`,
    async (failed) => {
      analysis_mock.models = []
      if (failed) analysis_mock.error = `Failed to load /assets/analysis.js: HTTP 503`
      mount(PhononsPage, { target: document.body })
      if (failed) {
        await vi.waitFor(() =>
          expect(doc_query(`.diagnostics-error details`).textContent).toContain(
            analysis_mock.error,
          ),
        )
        const reload = vi.spyOn(location, `reload`).mockImplementation(() => {})
        doc_query<HTMLButtonElement>(`.diagnostics-error button`).click()
        expect(reload).toHaveBeenCalledOnce()
        return
      }

      await vi.waitFor(() =>
        expect(document.body.textContent).toContain(
          `No per-material phonon diagnostics available for`,
        ),
      )
      expect(document.body.textContent).not.toContain(
        `Loading per-material phonon diagnostics`,
      )
      expect(document.querySelector(`.diagnostics-grid`)).toBeNull()
      expect(document.querySelector(`.diagnostics-error`)).toBeNull()
    },
  )

  it(`shows only phonon leaderboard columns and rows with phonon metrics`, () => {
    mount(PhononsPage, { target: document.body })

    const leaderboard = doc_query(`section.full-bleed`)
    const headers = [...leaderboard.querySelectorAll(`th`)].map(header_name)
    for (const column of [
      `Model`,
      `Links`,
      `κSRME`,
      `κSRE`,
      `κSRD`,
      `κfailed`,
      `Im(ω)`,
      `W1(ω)`,
    ]) {
      expect(headers).toContain(column)
    }
    for (const column of [`F1`, `DAF`, `Acc`, `Prec`]) {
      expect(headers).not.toContain(column)
    }

    const kappa_col_indices = [`κSRME`, `κSRE`, `κSRD`, `κfailed`, `Im(ω)`].map(
      (header) => headers.indexOf(header),
    )
    expect(kappa_col_indices.every((column_idx) => column_idx >= 0)).toBe(true)
    const rows = [...leaderboard.querySelectorAll<HTMLTableRowElement>(`tbody tr`)]
    expect(rows).toHaveLength(phonon_leaderboard_count)
    for (const row of rows) {
      const cells = [...row.querySelectorAll(`td`)]
      for (const column_idx of kappa_col_indices) {
        expect(cells[column_idx].textContent?.trim()).not.toBe(`n/a`)
      }
    }
    expect(sorted_header()?.textContent).toContain(`κSRME`)
    expect(sorted_header()?.getAttribute(`aria-sort`)).toBe(`ascending`)
  })

  it(`orders Compare models by κSRME and updates the model URL`, async () => {
    await mount_with_url(PhononsPage, `http://localhost/benchmarks/phonons`)

    // sort control offers alphabetical, κSRME, and date-added modes
    expect([...kappa_sort_select().options].map((option) => option.value)).toEqual([
      `kappa`,
      `name`,
      `date`,
    ])

    // default order is κSRME ascending (best phonon models first); option labels are
    // `<model_name> (<κSRME>)`, so strip the suffix to look up each model's metric
    const current_selection = kappa_selected_model() ?? ``
    const options = await kappa_model_options()
    const srme_values = options.map((option) => {
      const label = option.textContent?.trim() ?? ``
      const model_name = label.replace(/ \([^()]*\)$/, ``)
      const model = ACTIVE_MODELS.find((candidate) => candidate.model_name === model_name)
      return model?.metrics?.phonons?.kappa_103?.κ_SRME ?? Infinity
    })
    expect(srme_values.filter(Number.isFinite).length).toBeGreaterThan(1)
    expect(srme_values).toEqual(
      [...srme_values].toSorted((value_1, value_2) => value_1 - value_2),
    )

    const target_model = [`chgnet-0.3.0`, `mace-mp-0`]
      .map((model_key) => ACTIVE_MODELS.find((model) => model.model_key === model_key))
      .find((model) => model && !current_selection.includes(model.model_name))
    if (!target_model) throw new Error(`No alternate kappa model found`)
    const target_option = options.find((option) =>
      option.textContent?.includes(target_model.model_name),
    )
    if (!target_option)
      throw new Error(`No kappa option found for ${target_model.model_name}`)
    target_option.click()
    await tick()

    expect(kappa_selected_model()).toContain(target_model.model_name)
    expect(new URL(location.href).searchParams.get(`model`)).toBe(target_model.model_key)
  })

  it.each([`name`, `constructor`])(
    `restores URL state with model_sort=%s`,
    async (mode) => {
      const url = `http://localhost/benchmarks/phonons?model=mace-mp-0&model_sort=${mode}&x=F1&y=κ_SRE&sort=κ_SRE&dir=asc&openness=OSOD&heatmap=0`
      await mount_with_url(PhononsPage, url)

      expect(kappa_sort_select().value).toBe(mode === `name` ? `name` : `kappa`)
      const mace_name = ACTIVE_MODELS.find(
        (model) => model.model_key === `mace-mp-0`,
      )?.model_name
      if (!mace_name) throw new Error(`mace-mp-0 not found in ACTIVE_MODELS`)
      expect(kappa_selected_model()).toContain(mace_name)
      expect(heading_texts()).toContain(`Model Comparison: κSRE vs F1`)

      await expect_table_url_state(`κSRE`, `ascending`)
      expect(filter_summary_badge(`Openness`)).toContain(`(1/4)`)
    },
  )
})
