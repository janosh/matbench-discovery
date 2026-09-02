import {
  as_phonon_dos,
  dos_per_atom,
  build_kappa_parity_series,
  get_kappa_parity_point,
  has_kappa_parity_model,
  kappa_parity_manifest,
  kappa_model_asset,
  kappa_parity_asset_url,
  kappa_structure,
  load_kappa_parity_base,
  load_kappa_parity_model,
} from '$lib/parity/kappa-parity'
import type { KappaParityBase, KappaParityModel } from '$lib/parity/kappa-parity'
import { clear_asset_cache } from '$lib/asset-loader'
import type { AnyStructure } from 'matterviz/structure'
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest'
import { gzipped_json_response, request_url } from '../index'

beforeEach(clear_asset_cache)

afterEach(() => {
  vi.unstubAllGlobals()
  clear_asset_cache()
})

const dummy_structure = { sites: [], lattice: {} } as unknown as AnyStructure

const base: KappaParityBase = {
  material_ids: [`mp-1`, `mp-2`, `mp-3`],
  formulas: [`Si2`, `Ge2`, `Sn2`],
  kappa_dft: [10, 5, null],
  n_sites: [8, 4, 2],
  spacegroups: [225, 186, 1],
  structures: { 'mp-1': dummy_structure },
  dft_dos: { 'mp-1': { frequencies: [0, 1, 2], densities: [0, 1, 0] } },
}

const model: KappaParityModel = {
  model_key: `test-model`,
  kappa_ml: [8, null, 6],
  ml_dos: { 'mp-1': { frequencies: [0, 1, 2], densities: [0, 0.5, 0] } },
}

const first_model_key = Object.keys(kappa_parity_manifest.model_assets)[0]
if (!first_model_key) throw new Error(`kappa parity manifest has no model assets`)

function manifest_sized_base(overrides: Partial<KappaParityBase> = {}): KappaParityBase {
  const row_count = kappa_parity_manifest.row_count
  return {
    material_ids: Array<string>(row_count).fill(`mp-test`),
    formulas: Array<string>(row_count).fill(`Si2`),
    kappa_dft: Array<number | null>(row_count).fill(1),
    n_sites: Array<number | null>(row_count).fill(2),
    spacegroups: Array<number | null>(row_count).fill(225),
    structures: {},
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

  it(`rescales a peak-normalized DOS to 3 modes per atom for thermal properties`, () => {
    // triangle DOS on a unit grid: trapezoid integral 2 -> scale factor 3/2
    const per_atom = dos_per_atom({
      type: `phonon`,
      frequencies: [0, 1, 2, 3],
      densities: [0, 1, 1, 0],
    })
    expect(per_atom.frequencies).toEqual([0, 1, 2, 3])
    expect(per_atom.densities).toEqual([0, 1.5, 1.5, 0])
    const integral = per_atom.densities
      .slice(1)
      .reduce((sum, density, idx) => sum + (density + per_atom.densities[idx]) / 2, 0)
    expect(integral).toBeCloseTo(3, 12)
    expect(() =>
      dos_per_atom({ type: `phonon`, frequencies: [0, 1], densities: [0, 0] }),
    ).toThrow(`integrates to 0`)
  })

  it(`returns prebuilt structures and null for missing materials`, () => {
    expect(kappa_structure(base, `mp-1`)).toBe(dummy_structure)
    expect(kappa_structure(base, `mp-2`)).toBeNull()
  })

  it(`returns null (not throws) for an unparsable structure payload`, () => {
    // runs inside a $derived in the component, so a throw would crash the plot
    const bad_base = { ...base, structures: { 'mp-x': `not a valid structure` } }
    vi.spyOn(console, `error`).mockImplementation(() => {})
    vi.spyOn(console, `warn`).mockImplementation(() => {})
    expect(kappa_structure(bad_base, `mp-x`)).toBeNull()
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
    const base_url = kappa_parity_asset_url(kappa_parity_manifest.base.asset)
    const model_url = kappa_parity_asset_url(kappa_model_asset(first_model_key))
    const fetch_mock = vi.fn((url: RequestInfo | URL) => {
      const href = request_url(url)
      if (href === base_url) return gzipped_json_response(valid_base)
      if (href === model_url) return gzipped_json_response({ model: valid_model })
      return Promise.resolve(new Response(`missing`, { status: 404 }))
    })
    vi.stubGlobal(`fetch`, fetch_mock)

    await expect(load_kappa_parity_base()).resolves.toEqual(valid_base)
    await expect(load_kappa_parity_model(first_model_key)).resolves.toEqual(valid_model)
  })

  it.each([`kappa_dft`, `n_sites`, `spacegroups`] as const)(
    `rejects base %s with the wrong row count`,
    async (field) => {
      vi.stubGlobal(
        `fetch`,
        vi.fn(() => gzipped_json_response(manifest_sized_base({ [field]: [1] }))),
      )
      await expect(load_kappa_parity_base()).rejects.toThrow(
        `Invalid kappa parity ${field}: expected ${kappa_parity_manifest.row_count} rows`,
      )
    },
  )

  it(`rejects model data with the wrong row count`, async () => {
    vi.stubGlobal(
      `fetch`,
      vi.fn(() =>
        gzipped_json_response({ model: manifest_sized_model({ kappa_ml: [1] }) }),
      ),
    )
    await expect(load_kappa_parity_model(first_model_key)).rejects.toThrow(
      `Invalid kappa parity ${first_model_key}.kappa_ml: expected ${kappa_parity_manifest.row_count} rows`,
    )
  })
})
