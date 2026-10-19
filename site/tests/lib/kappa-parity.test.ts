import {
  as_phonon_dos,
  build_kappa_parity_series,
  get_kappa_parity_point,
  has_kappa_parity_model,
  kappa_parity_manifest,
  kappa_model_asset,
  kappa_parity_asset_url,
  load_kappa_parity_base,
  load_kappa_parity_model,
} from '#lib/parity/kappa-parity.js'
import type { KappaParityBase, KappaParityModel } from '#lib/parity/kappa-parity.js'
import * as kappa_parity from '#lib/parity/kappa-parity.js'
import * as kappa_modes from '#lib/parity/kappa-modes.js'
import KappaParityPlot from '#lib/plot/KappaParityPlot.svelte'
import kappa_103_analysis from '#figs/kappa-103-analysis.jsonl'
import { ACTIVE_MODELS, MODELS } from '#lib/models.svelte.js'
import { tick } from 'svelte'
import { describe, expect, it, onTestFinished, vi } from 'vite-plus/test'
import {
  doc_query,
  get_scatter_plot_props,
  gzipped_json_response,
  mount,
  stub_json_routes,
} from '../index.js'

const plot_mocks = vi.hoisted(() => ({
  ScatterPlot: vi.fn(),
  Dos: vi.fn(),
  PhononThermalPlot: vi.fn(),
  density_divisor: vi.fn(() => 2),
  spectral_loaded: vi.fn(),
}))
vi.mock(`matterviz/plot`, () => ({ ScatterPlot: plot_mocks.ScatterPlot }))
vi.mock(`matterviz/spectral`, () => {
  plot_mocks.spectral_loaded()
  return plot_mocks
})

const base: KappaParityBase = {
  material_ids: [`mp-1`, `mp-2`, `mp-3`],
  formulas: [`Si2`, `Ge2`, `Sn2`],
  kappa_dft: [10, 5, null],
  n_sites: [8, 4, 2],
  spacegroups: [225, 186, 1],
  dft_dos: { 'mp-1': { frequencies: [0, 1, 2], densities: [0, 1, 0] } },
}

const model: KappaParityModel = {
  model_key: `test-model`,
  kappa_ml: [8, null, 6],
  ml_dos: { 'mp-1': { frequencies: [0, 1, 2], densities: [0, 0.5, 0] } },
}

const first_model_key = Object.keys(kappa_parity_manifest.model_assets)[0]
if (!first_model_key) throw new Error(`kappa parity manifest has no model assets`)

it(`loads spectral plots only after selecting a material and renders both DOS sources`, async () => {
  vi.spyOn(kappa_modes, `has_kappa_modes`).mockReturnValue(true)
  const modes_load = vi
    .spyOn(kappa_modes, `load_kappa_modes`)
    .mockRejectedValueOnce(new Error(`HTTP 503: modes asset`))
    .mockResolvedValue({ model_key: model.model_key, materials: {} })
  vi.spyOn(kappa_parity, `load_kappa_parity_base`).mockResolvedValue(base)
  vi.spyOn(kappa_parity, `load_kappa_parity_model`)
    .mockRejectedValueOnce(new Error(`HTTP 503: kappa asset`))
    .mockResolvedValue(model)
  vi.spyOn(kappa_parity, `load_kappa_srme_map`).mockResolvedValue(
    new Map([[`mp-1`, 0.4]]),
  )
  mount(KappaParityPlot, {
    target: document.body,
    props: { model: { ...MODELS[0], model_key: first_model_key } },
  })
  await vi.waitFor(() =>
    expect(doc_query(`[role="alert"]`).textContent).toContain(`Could not load`),
  )
  expect(doc_query(`details`).textContent).toContain(`HTTP 503: kappa asset`)
  doc_query<HTMLButtonElement>(`.plot-state button`).click()
  const props = (await vi.waitFor(() =>
    get_scatter_plot_props(plot_mocks.ScatterPlot),
  )) as {
    series: { color_values: (number | null)[] }[]
    on_point_click: (event: { point: { series_idx: number; point_idx: number } }) => void
  }
  const color_values = () =>
    (get_scatter_plot_props(plot_mocks.ScatterPlot) as typeof props).series[0]
      .color_values
  expect(color_values()).toEqual([0.4])
  const color_select = doc_query<HTMLSelectElement>(`select[aria-label="Color metric"]`)
  color_select.value = `sre`
  // happy-dom lacks option:checked, which Svelte reads once in the change handler.
  vi.spyOn(color_select, `querySelector`).mockReturnValueOnce(
    color_select.selectedOptions[0],
  )
  color_select.dispatchEvent(new Event(`change`, { bubbles: true }))
  await tick()
  expect(color_values()).toEqual([(2 * Math.abs(8 - 10)) / (8 + 10)])
  expect(plot_mocks.spectral_loaded).not.toHaveBeenCalled()
  props.on_point_click({ point: { series_idx: 0, point_idx: 0 } })
  await vi.waitFor(() => expect(plot_mocks.Dos).toHaveBeenCalledOnce())
  expect(plot_mocks.PhononThermalPlot).toHaveBeenCalledTimes(2)
  // thermal plots get each DOS rescaled to 3 modes per atom over its integral (stubbed: 2)
  expect(
    plot_mocks.PhononThermalPlot.mock.calls.map(
      ([, thermal_props]) => thermal_props.dos.densities,
    ),
  ).toEqual([
    [0, 1.5, 0],
    [0, 0.75, 0],
  ])
  expect(document.querySelector(`.detail-panel`)?.textContent).toContain(`mp-1`)
  await vi.waitFor(() =>
    expect(doc_query(`.modes [role="alert"] details`).textContent).toContain(
      `HTTP 503: modes asset`,
    ),
  )
  doc_query<HTMLButtonElement>(`.modes .plot-state button`).click()
  await vi.waitFor(() =>
    expect(doc_query(`.modes`).textContent).toContain(
      `No phonon modes stored for this material`,
    ),
  )
  expect(modes_load).toHaveBeenCalledTimes(2)
  expect(document.querySelector(`.modes [role="alert"]`)).toBeNull()
  doc_query<HTMLButtonElement>(`button[aria-label="Close"]`).click()
  await tick()
  expect(document.querySelector(`.detail-panel`)).toBeNull()
})

it.each([`analysis`, `spectral`])(
  `offers a page reload after the %s module fails`,
  async (module) => {
    const failure = new Error(`Failed to load ${module}.js`)
    vi.spyOn(kappa_parity, `load_kappa_parity_base`).mockResolvedValue(base)
    vi.spyOn(kappa_parity, `load_kappa_parity_model`).mockResolvedValue(model)
    const srme_load = vi
      .spyOn(kappa_parity, `load_kappa_srme_map`)
      .mockResolvedValue(undefined)
    if (module === `analysis`) srme_load.mockRejectedValue(failure)
    else {
      vi.doMock(`matterviz/spectral`, () => {
        throw failure
      })
      onTestFinished(() => {
        vi.doMock(`matterviz/spectral`, () => plot_mocks)
      })
    }
    mount(KappaParityPlot, {
      target: document.body,
      props: { model: { ...MODELS[0], model_key: first_model_key } },
    })
    if (module === `spectral`) {
      const props = (await vi.waitFor(() =>
        get_scatter_plot_props(plot_mocks.ScatterPlot),
      )) as {
        on_point_click: (event: {
          point: { series_idx: number; point_idx: number }
        }) => void
      }
      props.on_point_click({ point: { series_idx: 0, point_idx: 0 } })
    }
    await vi.waitFor(() =>
      expect(doc_query(`[role="alert"]`).textContent).toContain(
        module === `analysis` ? failure.message : `Could not load thermal plots`,
      ),
    )
    // Vitest wraps deliberate module-factory failures; either path must expose details.
    expect(doc_query(`[role="alert"] details`).textContent).not.toBe(`Asset details`)
    const reload = vi.spyOn(location, `reload`).mockImplementation(() => {})
    const button = doc_query<HTMLButtonElement>(`[role="alert"] button`)
    expect(button.textContent).toBe(`Reload page`)
    button.click()
    expect(reload).toHaveBeenCalledOnce()
  },
)

function manifest_sized_base(overrides: Partial<KappaParityBase> = {}): KappaParityBase {
  const row_count = kappa_parity_manifest.row_count
  return {
    material_ids: Array<string>(row_count).fill(`mp-test`),
    formulas: Array<string>(row_count).fill(`Si2`),
    kappa_dft: Array<number | null>(row_count).fill(1),
    n_sites: Array<number | null>(row_count).fill(2),
    spacegroups: Array<number | null>(row_count).fill(225),
    dft_dos: {},
    ...overrides,
  }
}

const manifest_sized_model = (
  overrides: Partial<KappaParityModel> = {},
): KappaParityModel => ({
  model_key: first_model_key,
  kappa_ml: Array<number | null>(kappa_parity_manifest.row_count).fill(1),
  ml_dos: {},
  ...overrides,
})

describe(`kappa parity data helpers`, () => {
  it(`builds the scatter dropping rows missing either conductivity, carrying n_sites and spacegroup`, () => {
    const series = build_kappa_parity_series(base, model)
    expect(series.x).toEqual([10])
    expect(series.y).toEqual([8])
    expect(series.points.map((pt) => pt.material_id)).toEqual([`mp-1`])
    expect(series.points[0]).toMatchObject({ n_sites: 8, spacegroup: 225 })
  })

  it.each([
    [10, 8, 0, 0.2222],
    [4, 6, 0, 0.4],
  ])(
    `computes symmetric relative error for kappa %s vs %s`,
    (dft, ml, row_idx, expected_sre) => {
      const single = get_kappa_parity_point(
        { ...base, kappa_dft: [dft] },
        { ...model, kappa_ml: [ml] },
        row_idx,
      )
      expect(single?.sre).toBeCloseTo(expected_sre, 3)
    },
  )

  it.each([
    { label: `dft zero`, dft: 0, ml: 8 },
    { label: `dft negative`, dft: -1, ml: 8 },
    { label: `ml zero`, dft: 10, ml: 0 },
    { label: `ml negative`, dft: 10, ml: -1 },
  ])(`drops non-positive conductivities ($label), log-scale safe`, ({ dft, ml }) => {
    const point = get_kappa_parity_point(
      { ...base, kappa_dft: [dft] },
      { ...model, kappa_ml: [ml] },
      0,
    )
    expect(point).toBeNull()
  })

  it(`only reports models present in the manifest`, () => {
    expect(has_kappa_parity_model(first_model_key)).toBe(true)
    expect(has_kappa_parity_model(`missing-model`)).toBe(false)
    expect(has_kappa_parity_model(undefined)).toBe(false)
  })

  it(`maps phonon DOS to matterviz shape, or null when absent`, () => {
    expect(as_phonon_dos(base.dft_dos[`mp-1`])).toEqual({
      type: `phonon`,
      frequencies: [0, 1, 2],
      densities: [0, 1, 0],
    })
    expect(as_phonon_dos(base.dft_dos[`mp-2`])).toBeNull()
    expect(as_phonon_dos(model.ml_dos[`mp-1`])?.type).toBe(`phonon`)
    expect(as_phonon_dos(model.ml_dos[`mp-3`])).toBeNull()
  })

  it(`maps model keys to per-model release assets`, () => {
    const asset = kappa_model_asset(first_model_key)
    expect(asset).toContain(`-model-`)
    expect(kappa_parity_asset_url(asset)).toBe(`/kappa-parity/assets/${asset}`)
    expect(() => kappa_model_asset(`missing-model`)).toThrow(
      `No kappa parity model asset for missing-model`,
    )
  })

  it(`loads base and per-model assets through the manifest`, async () => {
    const valid_base = manifest_sized_base()
    const valid_model = manifest_sized_model()
    const fetch_mock = stub_json_routes({
      [kappa_parity_asset_url(kappa_parity_manifest.base.asset)]: valid_base,
      [kappa_parity_asset_url(kappa_model_asset(first_model_key))]: {
        model: valid_model,
      },
    })

    await expect(load_kappa_parity_base()).resolves.toEqual(valid_base)
    await expect(load_kappa_parity_model(first_model_key)).resolves.toEqual(valid_model)
    expect(fetch_mock).toHaveBeenCalledTimes(2)
  })

  it.each([
    ...([`kappa_dft`, `n_sites`, `spacegroups`] as const).map((field) => ({
      field,
      load: load_kappa_parity_base,
      asset: (rows?: number[]) => manifest_sized_base(rows && { [field]: rows }),
      valid: manifest_sized_base(),
    })),
    {
      field: `${first_model_key}.kappa_ml`,
      load: () => load_kappa_parity_model(first_model_key),
      asset: (rows?: number[]) => ({
        model: manifest_sized_model(rows && { kappa_ml: rows }),
      }),
      valid: manifest_sized_model(),
    },
  ])(
    `rejects $field with the wrong row count, then retries`,
    async ({ field, load, asset, valid }) => {
      const fetch_mock = vi
        .fn(() => gzipped_json_response(asset()))
        .mockImplementationOnce(() => gzipped_json_response(asset([1])))
      vi.stubGlobal(`fetch`, fetch_mock)
      await expect(load()).rejects.toThrow(
        `Invalid kappa parity ${field}: expected ${kappa_parity_manifest.row_count} rows, got 1`,
      )
      await expect(load()).resolves.toEqual(valid)
      expect(fetch_mock).toHaveBeenCalledTimes(2)
    },
  )

  it(`kappa-103 analysis payload matches the models' aggregate phonon metrics`, () => {
    const n_materials = kappa_103_analysis.material_ids.length
    for (const values of [
      kappa_103_analysis.formulas,
      kappa_103_analysis.spg_nums,
      kappa_103_analysis.kappa_dft,
    ]) {
      expect(values).toHaveLength(n_materials)
    }
    expect(
      kappa_103_analysis.spg_nums.every(
        (spg_num) => Number.isInteger(spg_num) && spg_num >= 1 && spg_num <= 230,
      ),
    ).toBe(true)

    for (const entry of kappa_103_analysis.models) {
      for (const values of [
        entry.kappa_ml,
        entry.srme,
        entry.srme_censored,
        entry.imag_modes,
        entry.broken_sym,
        entry.max_steps,
        entry.freq_w1,
      ]) {
        expect(values).toHaveLength(n_materials)
      }
      const metrics = ACTIVE_MODELS.find(
        (candidate) => candidate.model_key === entry.model_key,
      )?.metrics?.phonons?.kappa_103
      if (!metrics) throw new Error(`missing phonon metrics for ${entry.model_key}`)
      const payload_failure_rate =
        entry.srme_censored.filter((value) => value === true).length / n_materials
      expect(payload_failure_rate).toBeCloseTo(metrics.κ_failure_rate, 4)

      const payload_imaginary_rate =
        entry.imag_modes.filter((value) => value === true).length / n_materials
      expect(payload_imaginary_rate).toBeCloseTo(metrics.imaginary_mode_rate, 4)
      if (metrics.spectrum_w1 === null) {
        expect(entry.freq_w1_mean).toBeNull()
      } else {
        expect(entry.freq_w1_mean).toBeCloseTo(metrics.spectrum_w1, 4)
      }
    }
  })
})
