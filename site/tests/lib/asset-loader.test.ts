import {
  clear_asset_cache,
  load_json_asset,
  load_model_asset,
  parity_asset_resolver,
} from '#lib/asset-loader.js'
import {
  energy_parity_asset_url,
  energy_parity_manifest,
  has_energy_parity_model,
  model_asset as energy_model_asset,
} from '#lib/parity/energy-parity.js'
import { has_kappa_modes } from '#lib/parity/kappa-modes.js'
import {
  has_kappa_parity_model,
  kappa_model_asset,
  kappa_parity_asset_url,
  kappa_parity_manifest,
} from '#lib/parity/kappa-parity.js'
import { gzipSync } from 'node:zlib'
import { describe, expect, it, vi } from 'vite-plus/test'
import { gzipped_json_response } from '../index.js'

type AssetKind = `parity` | `modes`
type ModelAssets = Record<
  string,
  Record<string, { asset: string } | undefined> | undefined
>

const local_manifest = (model_assets: ModelAssets) => ({
  local_asset_base_url: `/local/assets/`,
  model_assets,
})

const kappa_modes_resolver = parity_asset_resolver(
  `kappa modes`,
  kappa_parity_manifest,
  undefined,
  `modes`,
)

describe(`parity_asset_resolver`, () => {
  // the module-level resolvers each pages and loaders import (kappa modes only exports
  // has_model, so its asset helpers come from the same manifest and asset kind)
  it.each([
    {
      kind: `energy`,
      asset_kind: `parity` as AssetKind,
      model_assets: energy_parity_manifest.model_assets as ModelAssets,
      base_url: `/energy-parity/assets`,
      has_model: has_energy_parity_model,
      model_asset: energy_model_asset,
      asset_url: energy_parity_asset_url,
    },
    {
      kind: `kappa`,
      asset_kind: `parity` as AssetKind,
      model_assets: kappa_parity_manifest.model_assets as ModelAssets,
      base_url: `/kappa-parity/assets`,
      has_model: has_kappa_parity_model,
      model_asset: kappa_model_asset,
      asset_url: kappa_parity_asset_url,
    },
    {
      kind: `kappa modes`,
      asset_kind: `modes` as AssetKind,
      model_assets: kappa_parity_manifest.model_assets as ModelAssets,
      base_url: `/kappa-parity/assets`,
      has_model: has_kappa_modes,
      model_asset: kappa_modes_resolver.model_asset,
      asset_url: kappa_modes_resolver.asset_url,
    },
  ])(
    `$kind resolver maps manifest models to $asset_kind asset URLs`,
    ({ kind, asset_kind, model_assets, base_url, has_model, model_asset, asset_url }) => {
      const entries = Object.entries(model_assets)
      expect(entries.some(([, kinds]) => kinds?.[asset_kind]?.asset)).toBe(true)
      for (const [model_key, kinds] of entries) {
        const asset = kinds?.[asset_kind]?.asset
        expect(has_model(model_key)).toBe(Boolean(asset))
        if (!asset) continue
        expect(model_asset(model_key)).toBe(asset)
        expect(asset).toContain(`-${asset_kind === `parity` ? `model` : `modes`}-`)
        expect(asset_url(asset)).toBe(`${base_url}/${asset}`)
      }
      expect(has_model(undefined)).toBe(false)
      expect(has_model(`missing-model`)).toBe(false)
      expect(() => model_asset(`missing-model`)).toThrow(
        `No ${kind} parity model asset for missing-model`,
      )
    },
  )

  it.each([
    { env_base_url: undefined, asset: `a.json.gz`, expected: `/local/assets/a.json.gz` },
    { env_base_url: ``, asset: `a.json.gz`, expected: `/local/assets/a.json.gz` },
    { env_base_url: `  `, asset: `/a.json.gz`, expected: `/local/assets/a.json.gz` },
    {
      env_base_url: ` https://cdn.test/assets/ `,
      asset: `a.json.gz`,
      expected: `https://cdn.test/assets/a.json.gz`,
    },
    {
      env_base_url: `https://cdn.test/assets//`,
      asset: `//nested/a.json.gz`,
      expected: `https://cdn.test/assets/nested/a.json.gz`,
    },
  ])(
    `env base URL $env_base_url joins $asset as $expected`,
    ({ env_base_url, asset, expected }) => {
      const { asset_url } = parity_asset_resolver(
        `test`,
        local_manifest({}),
        env_base_url,
      )
      expect(asset_url(asset)).toBe(expected)
    },
  )

  it.each([
    { asset_kind: `parity` as AssetKind, has_asset: true },
    { asset_kind: `modes` as AssetKind, has_asset: false },
  ])(
    `resolves only the requested $asset_kind asset kind`,
    ({ asset_kind, has_asset }) => {
      const resolver = parity_asset_resolver(
        `test`,
        local_manifest({
          'with-parity': { parity: { asset: `parity.json.gz` } },
          'empty-asset': { parity: { asset: `` }, modes: { asset: `` } },
        }),
        undefined,
        asset_kind,
      )
      expect(resolver.has_model(`with-parity`)).toBe(has_asset)
      // an empty asset name is as unusable as a missing one
      expect(resolver.has_model(`empty-asset`)).toBe(false)
      expect(() => resolver.model_asset(`empty-asset`)).toThrow(
        `No test parity model asset for empty-asset`,
      )
      if (has_asset) expect(resolver.model_asset(`with-parity`)).toBe(`parity.json.gz`)
      else expect(() => resolver.model_asset(`with-parity`)).toThrow(`with-parity`)
    },
  )
})

describe(`load_json_asset`, () => {
  it(`retries asset loads after a failed request`, async () => {
    const body = gzipSync(JSON.stringify({ ok: true }))
    const fetch_mock = vi
      .fn()
      .mockResolvedValueOnce(new Response(`missing`, { status: 503 }))
      .mockResolvedValueOnce(new Response(body, { status: 200 }))
    vi.stubGlobal(`fetch`, fetch_mock)

    const url = `/retry-${crypto.randomUUID()}.json.gz`
    await expect(load_json_asset(url)).rejects.toThrow(`HTTP 503`)
    await expect(load_json_asset(url)).resolves.toEqual({ ok: true })
    expect(fetch_mock).toHaveBeenCalledTimes(2)
  })

  it(`keeps a replacement cached request when a cleared request fails later`, async () => {
    const stale_response = Promise.withResolvers<Response>()
    const fetch_mock = vi
      .fn()
      .mockReturnValueOnce(stale_response.promise)
      .mockResolvedValueOnce(new Response(JSON.stringify({ fresh: true })))
    vi.stubGlobal(`fetch`, fetch_mock)
    const url = `/replaced-asset.json.gz`
    const stale_result = load_json_asset(url).catch((error: unknown) => error)
    clear_asset_cache()
    const replacement = load_json_asset(url)
    await expect(replacement).resolves.toEqual({ fresh: true })
    stale_response.reject(new Error(`old request failed`))
    await expect(stale_result).resolves.toEqual(new Error(`old request failed`))
    expect(load_json_asset(url)).toBe(replacement)
    expect(fetch_mock).toHaveBeenCalledTimes(2)
  })
})

describe(`load_model_asset`, () => {
  const model_key = `test-model`
  const valid = { model_key, values: [1, 2] }

  it.each([
    // asset stores the model object directly instead of under the { model } key
    {
      case: `missing { model } wrapper`,
      asset: valid,
      error: `No test model ${model_key} in its asset`,
    },
    {
      case: `model_key mismatch`,
      asset: { model: { ...valid, model_key: `other-model` } },
      error: `Invalid test model: expected ${model_key}, got other-model`,
    },
    {
      case: `failed caller validation`,
      asset: { model: { ...valid, values: [] } },
      error: `no values`,
    },
  ])(
    `rejects $case, then retries and caches the valid asset`,
    async ({ asset, error }) => {
      const fetch_mock = vi
        .fn(() => gzipped_json_response({ model: valid }))
        .mockImplementationOnce(() => gzipped_json_response(asset))
      vi.stubGlobal(`fetch`, fetch_mock)
      const load = () =>
        load_model_asset<typeof valid>(
          `test`,
          `/test-model.json.gz`,
          model_key,
          (model) => {
            if (model.values.length === 0) throw new Error(`no values`)
          },
        )
      await expect(load()).rejects.toThrow(error)
      await expect(load()).resolves.toEqual(valid)
      // the third load hits the cache of the validated retry
      await expect(load()).resolves.toEqual(valid)
      expect(fetch_mock).toHaveBeenCalledTimes(2)
    },
  )
})
