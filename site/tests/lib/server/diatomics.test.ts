import { fetch_diatomics_data, read_md_per_system } from '#lib/server/predictions.js'
import { mkdtemp, rm, writeFile } from 'node:fs/promises'
import { tmpdir } from 'node:os'
import { gzipSync } from 'node:zlib'
import { describe, expect, it, onTestFinished, vi } from 'vite-plus/test'

// skip retry backoff
vi.mock(`node:timers/promises`, () => ({ setTimeout: () => Promise.resolve() }))

describe(`diatomics server data loader`, () => {
  it(`does not fall back to remote data when local pred_file is corrupt`, async () => {
    const tmp_dir = await mkdtemp(`${tmpdir()}/diatomics-test-`)
    onTestFinished(() => rm(tmp_dir, { recursive: true }))
    await writeFile(`${tmp_dir}/diatomics.json.gz`, `not gzip`)
    const fetch_fn = vi.fn<typeof fetch>()

    await expect(
      fetch_diatomics_data(
        {
          pred_file: {
            name: `diatomics.json.gz`,
            url: `https://figshare.com/files/123`,
          },
        },
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
      fetch_diatomics_data(
        { pred_file: { name: `missing.json.gz`, url: `https://figshare.com/files/123` } },
        { fetch_fn, max_attempts: 3 },
      ),
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
        fetch_diatomics_data(
          {
            pred_file: { name: `missing.json.gz`, url: `https://figshare.com/files/123` },
          },
          { fetch_fn, max_attempts },
        ),
      ).rejects.toThrow(`max_attempts must be a positive integer`)
      expect(fetch_fn).not.toHaveBeenCalled()
    },
  )

  const pred_file = { name: `missing.json.gz`, url: `https://figshare.com/files/123` }
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
})
