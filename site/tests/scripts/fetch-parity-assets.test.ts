import { fetch_parity_assets } from '../../scripts/fetch-parity-assets'
import { createHash } from 'node:crypto'
import {
  mkdir,
  mkdtemp,
  readFile,
  readdir,
  rm,
  unlink,
  writeFile,
} from 'node:fs/promises'
import { tmpdir } from 'node:os'
import { setTimeout as sleep } from 'node:timers/promises'
import { afterEach, beforeEach, expect, it, vi } from 'vitest'

vi.mock(`node:timers/promises`, () => {
  const wait_mock = vi.fn(async () => {})
  return { setTimeout: wait_mock, default: { setTimeout: wait_mock } }
})

const bodies = { energy: `verified energy payload`, kappa: `verified kappa payload` }
const entry_for = (body: string, prefix = `energy-parity-v1`) => {
  const sha256 = createHash(`sha256`).update(body).digest(`hex`)
  return { asset: `${prefix}-${sha256.slice(0, 16)}.json.gz`, sha256 }
}
const energy_entry = entry_for(bodies.energy)
let site_dir: string
const asset_dir = (kind = `energy`) => `${site_dir}/static/${kind}-parity/assets`
const asset_path = (entry = energy_entry) => `${asset_dir()}/${entry.asset}`
const manifest_path = (kind = `energy`) =>
  `${site_dir}/src/lib/parity/${kind}-parity-manifest.json`
const fetch_mock = vi.fn<typeof fetch>()

beforeEach(async () => {
  site_dir = await mkdtemp(`${tmpdir()}/parity-download-`)
  await mkdir(`${site_dir}/src/lib/parity`, { recursive: true })
  for (const [kind, body] of Object.entries(bodies)) {
    const prefix = `${kind}-parity-v1`
    const entry = entry_for(body, prefix)
    await writeFile(
      manifest_path(kind),
      JSON.stringify({ asset_prefix: prefix, base: entry, model_assets: {} }),
    )
    await mkdir(asset_dir(kind), { recursive: true })
    await writeFile(
      `${asset_dir(kind)}/${entry.asset}`,
      kind === `energy` ? `corrupt` : body,
    )
  }
  fetch_mock.mockReset().mockImplementation(async () => new Response(bodies.energy))
  vi.stubGlobal(`fetch`, fetch_mock)
  vi.mocked(sleep).mockClear()
})

afterEach(async () => {
  vi.unstubAllGlobals()
  await rm(site_dir, { recursive: true, force: true })
})

it(`skips verified files and downloads model assets and structure bundles`, async () => {
  await writeFile(asset_path(), bodies.energy)
  const entries = [entry_for(`model`), entry_for(`structures`)]
  await writeFile(
    manifest_path(),
    JSON.stringify({
      asset_prefix: `energy-parity-v1`,
      base: energy_entry,
      model_assets: { model: { full_test_set: entries[0] } },
      structure_bundles: [entries[1]],
    }),
  )
  fetch_mock.mockImplementation(async (url) => {
    if (typeof url !== `string`) throw new Error(`Expected a URL string`)
    return new Response(url.endsWith(entries[0].asset) ? `model` : `structures`)
  })

  expect(await fetch_parity_assets(site_dir)).toBe(2)
  expect(fetch_mock).toHaveBeenCalledTimes(2)
  expect(await readFile(asset_path(entries[0]), `utf8`)).toBe(`model`)
  expect(await readFile(asset_path(entries[1]), `utf8`)).toBe(`structures`)
  expect(await fetch_parity_assets(site_dir)).toBe(0)
  expect(fetch_mock).toHaveBeenCalledTimes(2)
})

it(`retains a corrupt cache entry until its verified replacement is ready`, async () => {
  const response = Promise.withResolvers<Response>()
  fetch_mock.mockReturnValue(response.promise)
  const downloading = fetch_parity_assets(site_dir)
  await vi.waitFor(() => expect(fetch_mock).toHaveBeenCalledOnce())
  expect(await readFile(asset_path(), `utf8`)).toBe(`corrupt`)
  response.resolve(new Response(bodies.energy))
  expect(await downloading).toBe(1)
  expect(await readFile(asset_path(), `utf8`)).toBe(bodies.energy)
  expect(await readdir(asset_dir())).toEqual([energy_entry.asset])
})

it.each([`network`, `HTTP`, `body`])(
  `retries transient %s failures per file`,
  async (kind) => {
    if (kind === `network`)
      fetch_mock.mockRejectedValueOnce(new TypeError(`network down`))
    else if (kind === `HTTP`)
      fetch_mock.mockResolvedValueOnce(new Response(``, { status: 503 }))
    else {
      const response = new Response()
      vi.spyOn(response, `arrayBuffer`).mockRejectedValueOnce(
        new Error(`connection lost`),
      )
      fetch_mock.mockResolvedValueOnce(response)
    }
    expect(await fetch_parity_assets(site_dir)).toBe(1)
    expect(fetch_mock).toHaveBeenCalledTimes(2)
    expect(sleep).toHaveBeenCalledExactlyOnceWith(1000)
    expect(await readFile(asset_path(), `utf8`)).toBe(bodies.energy)
  },
)

it.each([404, 503])(
  `rejects HTTP %s without replacing an existing file`,
  async (status) => {
    fetch_mock.mockImplementation(async () => new Response(``, { status }))
    await expect(fetch_parity_assets(site_dir)).rejects.toThrow(`HTTP ${status}`)
    expect(fetch_mock).toHaveBeenCalledTimes(status === 404 ? 1 : 4)
    expect(await readFile(asset_path(), `utf8`)).toBe(`corrupt`)
  },
)

it(`rejects an invalid checksum without retrying or overwriting the cache`, async () => {
  fetch_mock.mockResolvedValue(new Response(`wrong bytes`))
  await expect(fetch_parity_assets(site_dir)).rejects.toThrow(
    `expected ${energy_entry.sha256}, received ${entry_for(`wrong bytes`).sha256}`,
  )
  expect(fetch_mock).toHaveBeenCalledOnce()
  expect(await readFile(asset_path(), `utf8`)).toBe(`corrupt`)
  expect(await readdir(asset_dir())).toEqual([energy_entry.asset])
})

it.each([
  { asset: `../${energy_entry.asset}` },
  { asset: energy_entry.asset.replace(`energy`, `kappa`) },
  { sha256: `not-a-checksum` },
])(`rejects invalid manifest entry %j before downloading`, async (invalid) => {
  await writeFile(
    manifest_path(),
    JSON.stringify({
      asset_prefix: `energy-parity-v1`,
      base: { ...energy_entry, ...invalid },
      model_assets: {},
    }),
  )
  await expect(fetch_parity_assets(site_dir)).rejects.toThrow(
    `Invalid energy parity asset`,
  )
  expect(fetch_mock).not.toHaveBeenCalled()
})

it(`propagates filesystem errors instead of treating them as missing cache files`, async () => {
  await unlink(asset_path())
  await mkdir(asset_path())
  await expect(fetch_parity_assets(site_dir)).rejects.toThrow(`EISDIR`)
  expect(fetch_mock).not.toHaveBeenCalled()
})
