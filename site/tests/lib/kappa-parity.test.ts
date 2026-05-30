import {
  build_kappa_parity_series,
  clear_kappa_parity_asset_cache,
  dft_phonon_dos,
  get_kappa_parity_point,
  has_kappa_parity_model,
  kappa_model_asset,
  kappa_parity_asset_url,
  kappa_structure,
  load_kappa_parity_base,
  load_kappa_parity_model,
  ml_phonon_dos,
} from '$lib/kappa-parity'
import type { KappaParityBase, KappaParityModel } from '$lib/kappa-parity'
import { kappa_parity_manifest } from '$lib/kappa-parity-manifest'
import type { AnyStructure } from 'matterviz/structure'
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest'
import { gzipped_json_response, request_url } from '../index'

beforeEach(clear_kappa_parity_asset_cache)

afterEach(() => {
  vi.unstubAllGlobals()
  clear_kappa_parity_asset_cache()
})

const dummy_structure = { sites: [], lattice: {} } as unknown as AnyStructure

const base: KappaParityBase = {
  material_ids: [`mp-1`, `mp-2`, `mp-3`],
  formulas: [`Si2`, `Ge2`, `Sn2`],
  kappa_dft: [10, 5, null],
  structures: { 'mp-1': dummy_structure },
  dft_dos: { 'mp-1': { frequencies: [0, 1, 2], densities: [0, 1, 0] } },
}

const model: KappaParityModel = {
  model_key: `test-model`,
  model_label: `Test Model`,
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
    structures: {},
    dft_dos: {},
    ...overrides,
  }
}

const manifest_sized_model = (
  overrides: Partial<KappaParityModel> = {},
): KappaParityModel => ({
  model_key: first_model_key,
  model_label: `Test Model`,
  kappa_ml: Array<number | null>(kappa_parity_manifest.row_count).fill(1),
  ml_dos: {},
  ...overrides,
})

describe(`kappa parity data helpers`, () => {
  it(`builds the scatter dropping rows missing either conductivity`, () => {
    const series = build_kappa_parity_series(base, model)
    expect(series.x).toEqual([10])
    expect(series.y).toEqual([8])
    expect(series.points.map((pt) => pt.material_id)).toEqual([`mp-1`])
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

  it(`carries n_sites/spacegroup, defaulting to null for older assets`, () => {
    const enriched = { ...base, n_sites: [8, 4, 2], spacegroups: [225, 186, 1] }
    const point = get_kappa_parity_point(enriched, model, 0)
    expect(point?.n_sites).toBe(8)
    expect(point?.spacegroup).toBe(225)

    // base without the optional arrays (older asset) must not throw and yields null
    const bare = get_kappa_parity_point(base, model, 0)
    expect(bare?.n_sites).toBeNull()
    expect(bare?.spacegroup).toBeNull()
  })

  it(`only reports models present in the manifest`, () => {
    expect(has_kappa_parity_model(first_model_key)).toBe(true)
    expect(has_kappa_parity_model(`missing-model`)).toBe(false)
    expect(has_kappa_parity_model(undefined)).toBe(false)
  })

  it(`maps phonon DOS to matterviz shape, or null when absent`, () => {
    expect(dft_phonon_dos(base, `mp-1`)).toEqual({
      type: `phonon`,
      frequencies: [0, 1, 2],
      densities: [0, 1, 0],
    })
    expect(dft_phonon_dos(base, `mp-2`)).toBeNull()
    expect(ml_phonon_dos(model, `mp-1`)?.type).toBe(`phonon`)
    expect(ml_phonon_dos(model, `mp-3`)).toBeNull()
  })

  it(`returns prebuilt structures and null for missing materials`, () => {
    expect(kappa_structure(base, `mp-1`)).toBe(dummy_structure)
    expect(kappa_structure(base, `mp-2`)).toBeNull()
  })

  it(`returns null (not throws) for an unparsable structure payload`, () => {
    // runs inside a $derived in the component, so a throw would crash the plot
    const bad_base = { ...base, structures: { 'mp-x': `not a valid structure` } }
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

  it.each([
    {
      kind: `base`,
      response: () => gzipped_json_response(manifest_sized_base({ kappa_dft: [1] })),
      load: () => load_kappa_parity_base(),
      error: `Invalid kappa parity kappa_dft: expected ${kappa_parity_manifest.row_count} rows`,
    },
    {
      kind: `model`,
      response: () =>
        gzipped_json_response({ model: manifest_sized_model({ kappa_ml: [1] }) }),
      load: () => load_kappa_parity_model(first_model_key),
      error: `Invalid kappa parity ${first_model_key}.kappa_ml: ` +
        `expected ${kappa_parity_manifest.row_count} rows`,
    },
  ])(`rejects $kind assets with the wrong row count`, async ({ response, load, error }) => {
    vi.stubGlobal(`fetch`, vi.fn(response))
    await expect(load()).rejects.toThrow(error)
  })
})
