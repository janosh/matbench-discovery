import { fetch_diatomics_data } from '$lib/server/predictions'
import { mkdtemp, rm, writeFile } from 'node:fs/promises'
import { tmpdir } from 'node:os'
import { describe, expect, it, onTestFinished, vi } from 'vitest'

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
})
