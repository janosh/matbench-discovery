import { clear_asset_cache } from '$lib/asset-loader'
import {
  build_phonon_mode_data,
  has_kappa_modes,
  load_kappa_modes,
} from '$lib/parity/kappa-modes'
import type { KappaModesMaterial, KappaModesModel } from '$lib/parity/kappa-modes'
import { kappa_parity_manifest } from '$lib/parity/kappa-parity'
import {
  nearest_qpoint_with_eigenvector,
  phonon_band_structure_from_modes,
} from 'matterviz/spectral'
import type { Complex, PhononQPointModes } from 'matterviz/spectral'
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest'
import { gzipped_json_response } from '../index'

beforeEach(clear_asset_cache)
afterEach(() => {
  vi.unstubAllGlobals()
  clear_asset_cache()
})

// two-atom cell, 6 bands, 3 q-points along GAMMA-X; eigenvectors stored only at the
// endpoints (index 0 and 2) like the site asset generator emits
const n_atoms = 2
const n_bands = 3 * n_atoms
const unit_eigenvector = (): Complex[][] =>
  Array.from({ length: n_atoms }, () =>
    Array.from({ length: 3 }, () => [1 / Math.sqrt(n_bands), 0]),
  )
const eigenvectors_at = () => Array.from({ length: n_bands }, unit_eigenvector)
const material = (): KappaModesMaterial => ({
  primitive: {
    lattice: [
      [4, 0, 0],
      [0, 4, 0],
      [0, 0, 4],
    ],
    symbols: [`Na`, `Cl`],
    masses: [22.99, 35.45],
    frac_coords: [
      [0, 0, 0],
      [0.5, 0.5, 0.5],
    ],
  },
  segments: [{ start_index: 0, end_index: 2, start_label: `GAMMA`, end_label: `X` }],
  q_points: [
    [0, 0, 0],
    [0.25, 0, 0],
    [0.5, 0, 0],
  ],
  distances: [0, 0.5, 1],
  frequencies: Array.from({ length: 3 }, (_, q_idx) =>
    Array.from({ length: n_bands }, (_unused, band_idx) => q_idx + band_idx / 10),
  ),
  eigenvectors: { 0: eigenvectors_at(), 2: eigenvectors_at() },
})

describe(`build_phonon_mode_data`, () => {
  it(`converts the compact asset record into matterviz PhononModeData`, () => {
    const data = build_phonon_mode_data(material())
    expect(data.n_atoms).toBe(2)
    expect(data.atoms.map((atom) => atom.symbol)).toEqual([`Na`, `Cl`])
    expect(data.atoms[1]).toEqual({
      symbol: `Cl`,
      mass: 35.45,
      coordinates: [0.5, 0.5, 0.5],
    })
    expect(data.lattice).toEqual([
      [4, 0, 0],
      [0, 4, 0],
      [0, 0, 4],
    ])
    expect(data.path_segments).toEqual([
      { start_index: 0, end_index: 2, start_label: `GAMMA`, end_label: `X` },
    ])
    expect(data.qpoints).toHaveLength(3)
    expect(data.qpoints[1]).toMatchObject({ q_position: [0.25, 0, 0], distance: 0.5 })
    expect(data.qpoints[1].modes.map((mode) => mode.frequency)).toEqual([
      1, 1.1, 1.2, 1.3, 1.4, 1.5,
    ])
    // eigenvectors only at the labeled endpoints, null in between
    const has_eigenvectors = (qpoint: PhononQPointModes) =>
      qpoint.modes.every((mode) => mode.eigenvector !== null)
    expect(data.qpoints.map(has_eigenvectors)).toEqual([true, false, true])
    expect(data.qpoints[0].modes[0].eigenvector).toEqual(unit_eigenvector())
    expect(nearest_qpoint_with_eigenvector(data, 1, 0)).toBeOneOf([0, 2])
    // the explorer's band view derives from the same data, so it must accept it
    const bands = phonon_band_structure_from_modes(data)
    expect(bands.qpoints).toHaveLength(3)
  })
})

describe(`kappa modes assets`, () => {
  it(`reports mode assets only for models listed in the manifest`, () => {
    expect(has_kappa_modes(undefined)).toBe(false)
    expect(has_kappa_modes(`missing-model`)).toBe(false)
    for (const model_key of Object.keys(kappa_parity_manifest.mode_assets)) {
      expect(has_kappa_modes(model_key)).toBe(true)
      // every model with modes also has the parity asset the page loads first
      expect(model_key in kappa_parity_manifest.model_assets).toBe(true)
    }
  })

  it(`throws for models without a mode asset instead of fetching`, async () => {
    const fetch_mock = vi.fn()
    vi.stubGlobal(`fetch`, fetch_mock)
    await expect(load_kappa_modes(`missing-model`)).rejects.toThrow(
      `No kappa modes parity model asset for missing-model`,
    )
    expect(fetch_mock).not.toHaveBeenCalled()
  })

  it(`loads and validates a mode asset from the manifest`, async () => {
    // no shipped model carries modes yet, so register one in the (mutable) manifest
    // object; the resolver reads it at call time
    const model_key = `modes-test-model`
    const mode_assets: Record<string, { asset: string; sha256: string } | undefined> =
      kappa_parity_manifest.mode_assets
    mode_assets[model_key] = {
      asset: `${kappa_parity_manifest.asset_prefix}-modes-${model_key}-${`0`.repeat(16)}.json.gz`,
      sha256: `a`.repeat(64),
    }
    try {
      expect(has_kappa_modes(model_key)).toBe(true)
      const payload: KappaModesModel = {
        model_key,
        materials: { 'mp-1': material() },
      }
      const fetch_mock = vi.fn(() => gzipped_json_response({ model: payload }))
      vi.stubGlobal(`fetch`, fetch_mock)
      await expect(load_kappa_modes(model_key)).resolves.toEqual(payload)
      expect(fetch_mock).toHaveBeenCalledWith(
        `/kappa-parity/assets/${mode_assets[model_key]?.asset}`,
      )
      clear_asset_cache()
      vi.stubGlobal(
        `fetch`,
        vi.fn(() =>
          gzipped_json_response({ model: { ...payload, model_key: `other-model` } }),
        ),
      )
      await expect(load_kappa_modes(model_key)).rejects.toThrow(
        `expected ${model_key}, got other-model`,
      )
    } finally {
      delete mode_assets[model_key]
    }
  })
})
