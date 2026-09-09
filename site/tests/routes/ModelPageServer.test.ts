import type { ModelData } from '$lib/types'
import { entries, load } from '$routes/models/[slug]/+page.server'
import { mkdtemp, rm, writeFile } from 'node:fs/promises'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { gzipSync } from 'node:zlib'
import { afterEach, describe, expect, it, onTestFinished, vi } from 'vitest'

const model: Pick<ModelData, `model_key` | `metrics`> = vi.hoisted(() => ({
  model_key: `model-key`,
}))

vi.mock(`$lib/models.svelte`, () => ({ MODELS: [model] }))

const load_model = () =>
  load({ params: { slug: model.model_key } } as Parameters<typeof load>[0])

// CSVs written by scripts/evals/md.py use eV and empty cells for unavailable values.
const csv_fixture = [
  `\uFEFFsystem,energy_rmse,force_rmse,pressure_error,provenance,note`,
  `CsSnI3_500K_Ivor_VASP,0.0021,0.04175,,"{""code"":""abc,123""}","two ""quoted"" lines\nwith a comma, too"`,
  `bulkCu_1000K,,,45.6,,`,
].join(`\r\n`)

describe(`model-key routes`, () => {
  afterEach(() => {
    vi.unstubAllGlobals()
  })

  it.each([undefined, `missing.csv.gz`])(
    `resolves the canonical key with unavailable MD predictions (%s)`,
    async (name) => {
      model.metrics = { md: { pred_file: name ? { name } : undefined } }
      expect(entries()).toEqual([{ slug: `model-key` }])
      expect(await load_model()).toEqual({ model, md_per_system: null })
    },
  )

  it.each([
    [false, false],
    [false, true],
    [true, false],
    [true, true],
  ])(
    `loads quoted MD rows with RMSEs in meV (remote=%s, gzip=%s)`,
    async (remote, gzip) => {
      const tmp_dir = await mkdtemp(join(tmpdir(), `md-per-system-`))
      onTestFinished(() => rm(tmp_dir, { recursive: true }))
      const csv_path = join(tmp_dir, gzip ? `metrics.csv.gz` : `metrics.csv`)
      const bytes = gzip ? gzipSync(csv_fixture) : csv_fixture
      if (!remote) await writeFile(csv_path, bytes)
      const fetch_mock = vi.fn<typeof fetch>().mockResolvedValue(new Response(bytes))
      vi.stubGlobal(`fetch`, fetch_mock)
      model.metrics = {
        md: { pred_file: { name: csv_path, url: `https://figshare.com/files/123` } },
      }

      expect(await load_model()).toEqual({
        model,
        md_per_system: [
          {
            system: `CsSnI3_500K_Ivor_VASP`,
            provenance: `{"code":"abc,123"}`,
            energy_rmse: 2.1,
            force_rmse: 41.75,
            note: `two "quoted" lines\nwith a comma, too`,
          },
          {
            system: `bulkCu_1000K`,
            pressure_error: 45.6,
          },
        ],
      })
      expect(fetch_mock).toHaveBeenCalledTimes(remote ? 1 : 0)
      if (remote) {
        expect(fetch_mock).toHaveBeenCalledWith(
          `https://ndownloader.figshare.com/files/123`,
          { signal: expect.any(AbortSignal) },
        )
      }
    },
  )

  it.each([
    [202, `challenge`],
    [404, `missing`],
    [500, `failed`],
    [200, ``],
    [200, `system,value\nCu,"unterminated`],
    [200, `system,value\nCu,"quoted"trailing`],
    [200, `system,value\nCu,1,2`],
    [200, `system,value\nCu`],
    [200, `system,system\nCu,1`],
  ])(`rejects unusable published predictions (%s, %s)`, async (status, body) => {
    const fetch_mock = vi
      .fn<typeof fetch>()
      .mockResolvedValue(new Response(body, { status }))
    vi.stubGlobal(`fetch`, fetch_mock)
    model.metrics = {
      md: { pred_file: { name: `missing.csv.gz`, url: `https://example.org/md.csv.gz` } },
    }
    await expect(load_model()).rejects.toThrow(
      status === 200 ? `Invalid MD prediction CSV` : `HTTP ${status}`,
    )
    expect(fetch_mock).toHaveBeenCalledTimes(1)
  })

  it(`rejects unknown keys`, async () => {
    await expect(
      load({
        params: { slug: `missing` },
      } as Parameters<typeof load>[0]),
    ).rejects.toMatchObject({ status: 404 })
  })
})
