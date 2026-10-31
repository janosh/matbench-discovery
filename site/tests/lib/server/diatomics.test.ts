import { by_benchmark_added_desc } from '#lib'
import { ACTIVE_MODELS, has_diatomics_curves } from '#lib/models.svelte.js'
import * as predictions from '#lib/server/predictions.js'
import { fetch_diatomics_data, read_md_per_system } from '#lib/server/predictions.js'
import { load } from '#routes/benchmarks/diatomics/+page.server.js'
import { createHash } from 'node:crypto'
import { mkdtemp, rm, writeFile } from 'node:fs/promises'
import { tmpdir } from 'node:os'
import { gzipSync } from 'node:zlib'
import { describe, expect, it, onTestFinished, vi } from 'vite-plus/test'

// skip retry backoff
vi.mock(`node:timers/promises`, () => ({ setTimeout: () => Promise.resolve() }))

const pred_file = { name: `missing.json.gz`, url: `https://figshare.com/files/123` }

describe(`diatomics server data loader`, () => {
  it(`does not fall back to remote data when local pred_file is corrupt`, async () => {
    const tmp_dir = await mkdtemp(`${tmpdir()}/diatomics-test-`)
    onTestFinished(() => rm(tmp_dir, { recursive: true }))
    await writeFile(`${tmp_dir}/diatomics.json.gz`, `not gzip`)
    const fetch_fn = vi.fn<typeof fetch>()

    await expect(
      fetch_diatomics_data(
        { pred_file: { ...pred_file, name: `diatomics.json.gz` } },
        { fetch_fn, root_dir: tmp_dir },
      ),
    ).rejects.toThrow(/did not return gzip-compressed JSON/)
    expect(fetch_fn).not.toHaveBeenCalled()
  })

  it(`does not retry Figshare WAF challenge responses`, async () => {
    const fetch_fn = vi.fn<typeof fetch>().mockResolvedValue(
      new Response(null, {
        status: 202,
        headers: {
          'content-type': `text/html; charset=UTF-8`,
          'x-amzn-waf-action': `challenge`,
        },
      }),
    )

    await expect(
      fetch_diatomics_data({ pred_file }, { fetch_fn, max_attempts: 3 }),
    ).rejects.toThrow(/Figshare WAF challenge: challenge/)
    expect(fetch_fn).toHaveBeenCalledTimes(1)
    // the WAF-gated figshare.com/files/<id> URL is rewritten to the ndownloader host
    expect(fetch_fn).toHaveBeenCalledWith(
      `https://ndownloader.figshare.com/files/123`,
      expect.objectContaining({ signal: expect.any(AbortSignal) }),
    )
  })

  it.each([0, -1, 1.5, Number.NaN])(
    `rejects invalid max_attempts value %s`,
    async (max_attempts) => {
      const fetch_fn = vi.fn<typeof fetch>()

      await expect(
        fetch_diatomics_data({ pred_file }, { fetch_fn, max_attempts }),
      ).rejects.toThrow(`max_attempts must be a positive integer`)
      expect(fetch_fn).not.toHaveBeenCalled()
    },
  )

  const loaders = {
    diatomics: (fetch_fn: typeof fetch) =>
      fetch_diatomics_data({ pred_file }, { fetch_fn, max_attempts: 3 }),
    md: (fetch_fn: typeof fetch) =>
      read_md_per_system(
        { metrics: { md: { pred_file } } },
        { fetch_fn, max_attempts: 3 },
      ),
  }
  const body = { diatomics: gzipSync(`{}`), md: `system,rdf_error\nCu,1.5` }

  it.each([
    [`diatomics`, [503, 500, 200], 3, {}],
    [`md`, [502, 200], 2, [{ system: `Cu`, rdf_error: 1.5 }]],
    [`diatomics`, [500, 503, 502], 3, /failed after 3 attempts.*HTTP 502/],
    [`md`, [404], 1, /HTTP 404/],
  ] as const)(
    `%s retries 5xx but not 4xx (statuses %j)`,
    async (loader, statuses, n_calls, outcome) => {
      const fetch_fn = vi.fn<typeof fetch>()
      for (const status of statuses) {
        fetch_fn.mockResolvedValueOnce(
          new Response(status === 200 ? body[loader] : null, { status }),
        )
      }
      const result = loaders[loader](fetch_fn)
      if (outcome instanceof RegExp) await expect(result).rejects.toThrow(outcome)
      else await expect(result).resolves.toEqual(outcome)
      expect(fetch_fn).toHaveBeenCalledTimes(n_calls)
    },
  )

  // local builds must never silently use data differing from the registered artifact
  const registered = `system\nCu`
  const figshare_url = `https://figshare.com/files/1`
  it.each([
    [registered, figshare_url, registered, 0, undefined],
    [`system\nAg`, figshare_url, registered, 1, undefined], // same size: only md5 catches it
    [`system\nAg`, undefined, registered, 0, `and no url is registered`],
    [
      `system\nAg`,
      figshare_url,
      `system\nAg`,
      1,
      `${figshare_url} differs from its registered artifact (md5`,
    ],
    [`system\nCu,extra`, undefined, registered, 0, `size 15 != registered 9`],
  ])(
    `pred_file: local %j, url %s, remote %j`,
    async (local, url, remote, n_fetches, error) => {
      const tmp_dir = await mkdtemp(`${tmpdir()}/pred-file-test-`)
      onTestFinished(() => rm(tmp_dir, { recursive: true }))
      await writeFile(`${tmp_dir}/pred-file`, local)
      const md5 = createHash(`md5`).update(registered).digest(`hex`)
      const file = { name: `pred-file`, size: registered.length, md5, url }
      const fetch_fn = vi.fn<typeof fetch>(async () => new Response(remote))
      const result = read_md_per_system(
        { metrics: { md: { pred_file: file } } },
        { fetch_fn, root_dir: tmp_dir },
      )
      if (error) await expect(result).rejects.toThrow(error)
      else await expect(result).resolves.toEqual([{ system: `Cu` }])
      expect(fetch_fn).toHaveBeenCalledTimes(n_fetches)
    },
  )

  const load_page = () => load({} as Parameters<typeof load>[0])
  const diatomics_keys = ACTIVE_MODELS.filter(has_diatomics_curves).map(
    ({ model_key }) => model_key,
  )

  it(`page load returns model keys (not ModelData) with homonuclear energies only`, async () => {
    const curves = {
      distances: [1, 2],
      'homo-nuclear': { 'H-H': { energies: [0, 1], forces: [[[0, 0, 0]]] } },
      'hetero-nuclear': { 'H-He': { energies: [0, 1] } },
    }
    vi.spyOn(predictions, `fetch_diatomics_data`).mockResolvedValue(
      curves as unknown as Awaited<ReturnType<typeof fetch_diatomics_data>>,
    )
    const data = await load_page()
    if (!data) throw new Error(`diatomics load returned nothing`)
    // newest first: the page's color palette and curve picker follow this order
    expect(data.diatomic_model_keys).toEqual(
      ACTIVE_MODELS.filter(has_diatomics_curves)
        .toSorted(by_benchmark_added_desc)
        .map(({ model_key }) => model_key),
    )
    expect(Object.keys(data).toSorted()).toEqual([
      `diatomic_curves`,
      `diatomic_model_keys`,
      `reference_names`,
    ])
    expect(data.diatomic_curves[diatomics_keys[0]]).toEqual({
      distances: [1, 2],
      'homo-nuclear': { 'H-H': { energies: [0, 1] } },
    })
  })

  it(`page load fails (and so the build) when any model's curves fail to download`, async () => {
    const failing_key = diatomics_keys[1]
    vi.spyOn(predictions, `fetch_diatomics_data`).mockImplementation(async (opts) => {
      if (
        opts.pred_file ===
        ACTIVE_MODELS.find(({ model_key }) => model_key === failing_key)?.metrics
          ?.diatomics?.pred_file
      )
        throw new Error(`HTTP 503`)
      return { distances: [], 'homo-nuclear': {} }
    })
    const result = load_page()
    await expect(result).rejects.toThrow(
      `Failed to fetch diatomics curves for ${failing_key}`,
    )
    await expect(result).rejects.toHaveProperty(`cause.message`, `HTTP 503`)
  })
})
