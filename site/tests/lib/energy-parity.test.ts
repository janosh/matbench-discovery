import {
  build_energy_parity_series,
  clear_energy_parity_asset_cache,
  energy_parity_asset_url,
  get_energy_parity_point,
  has_energy_parity_model,
  load_energy_parity_base,
  load_energy_parity_model,
  load_json_asset,
  load_wbm_structure,
  model_asset,
  structure_bundle_for_shard,
  structure_shard_idx,
  structure_popup_placement,
} from '$lib/energy-parity'
import type { EnergyParityBase, EnergyParityModel } from '$lib/energy-parity'
import { energy_parity_manifest } from '$lib/energy-parity-manifest'
import { gzipSync } from 'node:zlib'
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest'
import { gzipped_json_response, request_url } from '../index'

beforeEach(clear_energy_parity_asset_cache)

afterEach(() => {
  vi.unstubAllGlobals()
  clear_energy_parity_asset_cache()
})

const base: EnergyParityBase = {
  material_ids: [`wbm-1-1`, `wbm-1-2`, `wbm-1-3`],
  formulas: [`Ac6 U2`, `Ac1 Th3`, `Ac1 U3`],
  n_sites: [8, 4, 4],
  e_form_true: [0.5, 0.1, 0.7],
  each_true: [0.6, 0.2, 0.8],
  structure_shard_size: 2,
}

const model: EnergyParityModel = {
  model_key: `test-model`,
  model_label: `Test Model`,
  e_form_pred: [0.4, null, 0.9],
}

const first_model_key = Object.keys(energy_parity_manifest.model_assets)[0]
if (!first_model_key) throw new Error(`energy parity manifest has no model assets`)

const first_structure_bundle = energy_parity_manifest.structure_bundles[0]
if (!first_structure_bundle) {
  throw new Error(`energy parity manifest has no structure bundles`)
}

function manifest_sized_base(overrides: Partial<EnergyParityBase> = {}): EnergyParityBase {
  const row_count = energy_parity_manifest.row_count
  return {
    material_ids: Array<string>(row_count).fill(`wbm-1-test`),
    formulas: Array<string>(row_count).fill(`Te2`),
    n_sites: Array<number | null>(row_count).fill(2),
    e_form_true: Array<number | null>(row_count).fill(0),
    each_true: Array<number | null>(row_count).fill(0),
    structure_shard_size: energy_parity_manifest.structure_shard_size,
    ...overrides,
  }
}

const manifest_sized_model = (
  model_key = first_model_key,
  overrides: Partial<EnergyParityModel> = {},
): EnergyParityModel => ({
  model_key,
  model_label: `Test Model`,
  e_form_pred: Array<number | null>(energy_parity_manifest.row_count).fill(0),
  ...overrides,
})

describe(`energy parity data helpers`, () => {
  it(`builds typed scatter arrays and preserves source row ids`, () => {
    const series = build_energy_parity_series(base, model, `e-form`)

    expect([...series.x]).toEqual([0.5, expect.closeTo(0.7)])
    expect([...series.y]).toEqual([expect.closeTo(0.4), expect.closeTo(0.9)])
    expect([...series.point_ids]).toEqual([0, 2])
    expect([...series.size_values]).toEqual([8, 4])
  })

  it(`converts formation-energy predictions to convex-hull-distance predictions`, () => {
    const point = get_energy_parity_point(base, model, 2, `each`)

    expect(point).toMatchObject({
      row_idx: 2,
      material_id: `wbm-1-3`,
      formula: `Ac1 U3`,
      n_sites: 4,
      x: 0.8,
    })
    expect(point?.y).toBeCloseTo(1.0)
    expect(point?.error).toBeCloseTo(0.2)
  })

  it(`does not treat a stale injected model as loaded`, () => {
    expect(has_energy_parity_model(model, `test-model`)).toBe(true)
    expect(has_energy_parity_model(model, `other-model`)).toBe(false)
    expect(has_energy_parity_model(undefined, `test-model`)).toBe(false)
  })

  it(`loads gzipped JSON assets and caches duplicate requests`, async () => {
    const body = gzipSync(JSON.stringify({ ok: true }))
    const fetch_mock = vi.fn(() =>
      Promise.resolve(new Response(body, { status: 200 })),
    )
    vi.stubGlobal(`fetch`, fetch_mock)

    const url = `/fixture-${crypto.randomUUID()}.json.gz`
    await expect(load_json_asset(url)).resolves.toEqual({ ok: true })
    await expect(load_json_asset(url)).resolves.toEqual({ ok: true })
    expect(fetch_mock).toHaveBeenCalledTimes(1)
  })

  it(`retries asset loads after a failed request`, async () => {
    const body = gzipSync(JSON.stringify({ ok: true }))
    const fetch_mock = vi.fn()
      .mockResolvedValueOnce(new Response(`missing`, { status: 503 }))
      .mockResolvedValueOnce(new Response(body, { status: 200 }))
    vi.stubGlobal(`fetch`, fetch_mock)

    const url = `/retry-${crypto.randomUUID()}.json.gz`
    await expect(load_json_asset(url)).rejects.toThrow(`HTTP 503`)
    await expect(load_json_asset(url)).resolves.toEqual({ ok: true })
    expect(fetch_mock).toHaveBeenCalledTimes(2)
  })

  it(`maps model keys to readable per-model release assets`, () => {
    const asset = model_asset(first_model_key)
    const model_assets = energy_parity_manifest.model_assets as Record<
      string,
      { asset: string } | undefined
    >

    expect(asset).toBe(model_assets[first_model_key]?.asset)
    expect(asset).toContain(`-model-`)
    expect(energy_parity_asset_url(asset)).toBe(`/energy-parity/assets/${asset}`)
    expect(() => model_asset(`missing-model`)).toThrow(
      `No energy parity model asset for missing-model`,
    )
  })

  it(`loads base and per-model assets through the manifest`, async () => {
    const valid_base = manifest_sized_base()
    const valid_model = manifest_sized_model()
    const base_url = energy_parity_asset_url(energy_parity_manifest.base.asset)
    const model_url = energy_parity_asset_url(model_asset(first_model_key))
    const fetch_mock = vi.fn((url: RequestInfo | URL) => {
      const href = request_url(url)
      if (href === base_url) return gzipped_json_response(valid_base)
      if (href === model_url) return gzipped_json_response({ model: valid_model })
      return Promise.resolve(new Response(`missing`, { status: 404 }))
    })
    vi.stubGlobal(`fetch`, fetch_mock)

    await expect(load_energy_parity_base()).resolves.toEqual(valid_base)
    await expect(load_energy_parity_model(first_model_key)).resolves.toEqual(valid_model)
    expect(fetch_mock).toHaveBeenCalledTimes(2)
  })

  it.each([
    {
      kind: `base`,
      response: () => gzipped_json_response(manifest_sized_base({ formulas: [] })),
      load: () => load_energy_parity_base(),
      error: `Invalid energy parity formulas: expected ${energy_parity_manifest.row_count} rows`,
    },
    {
      kind: `model`,
      response: () =>
        gzipped_json_response({
          model: manifest_sized_model(first_model_key, { e_form_pred: [0] }),
        }),
      load: () => load_energy_parity_model(first_model_key),
      error: `Invalid energy parity ${first_model_key}.e_form_pred: ` +
        `expected ${energy_parity_manifest.row_count} rows`,
    },
  ])(`rejects $kind assets with the wrong row count`, async ({ response, load, error }) => {
    vi.stubGlobal(`fetch`, vi.fn(response))
    await expect(load()).rejects.toThrow(error)
  })

  it(`maps row indices to bundled structure shard assets`, () => {
    expect(structure_shard_idx(base, 0)).toBe(0)
    expect(structure_shard_idx(base, 1)).toBe(0) // mid-shard row (distinguishes floor from ceil)
    expect(structure_shard_idx(base, 2)).toBe(1)
    expect(structure_bundle_for_shard(1)).toBe(first_structure_bundle)
  })

  it(`loads structures from the shard nested inside a bundle`, async () => {
    const structure = { sites: [], lattice: [] }
    const material_id = base.material_ids[2] ?? ``
    const structure_url = energy_parity_asset_url(first_structure_bundle.asset)
    const fetch_mock = vi.fn((url: RequestInfo | URL) => {
      const href = request_url(url)
      if (href === structure_url) {
        return gzipped_json_response({
          shards: { 1: { [material_id]: structure } },
        })
      }
      return Promise.resolve(new Response(`missing`, { status: 404 }))
    })
    vi.stubGlobal(`fetch`, fetch_mock)

    await expect(load_wbm_structure(base, 2, material_id)).resolves.toEqual(structure)
    expect(fetch_mock).toHaveBeenCalledTimes(1)
  })

  it(`places structure popup in side gutters before falling back over the plot`, () => {
    const placement = {
      plot_width: 800,
      plot_height: 520,
      popup_width: 500,
    }

    expect(
      structure_popup_placement({
        ...placement,
        viewport_width: 1_600,
        plot_left: 700,
      }),
    ).toEqual({ side: `left`, place_right: false, left: 64, top: 260 })

    expect(
      structure_popup_placement({
        ...placement,
        viewport_width: 1_600,
        plot_left: 16,
      }),
    ).toEqual({ side: `right`, place_right: true, left: 776, top: 260 })

    expect(
      structure_popup_placement({
        ...placement,
        viewport_width: 900,
        plot_left: 100,
      }),
    ).toEqual({ side: `overlay`, place_right: false, left: 776, top: 260 })
  })
})
