import { createHash } from 'node:crypto'
import { mkdir, readFile, rename, writeFile } from 'node:fs/promises'
import { setTimeout as sleep } from 'node:timers/promises'
import { fileURLToPath } from 'node:url'

type Asset = { asset: string; sha256: string }
type Manifest = {
  asset_prefix: string
  base: Asset
  model_assets: Record<string, Record<string, Asset>>
  structure_bundles?: Asset[]
}

const sha256 = (bytes: Uint8Array): string =>
  createHash(`sha256`).update(bytes).digest(`hex`)

// Retry each transient failure separately, keeping successful downloads cached.
async function download(url: string): Promise<Uint8Array> {
  for (let attempt = 0; ; attempt++) {
    let response: Response | undefined
    try {
      response = await fetch(url, { signal: AbortSignal.timeout(120_000) })
      if (!response.ok) throw new Error(`HTTP ${response.status}`)
      return new Uint8Array(await response.arrayBuffer())
    } catch (error) {
      if (
        attempt === 3 ||
        (response &&
          !response.ok &&
          ![408, 429, 500, 502, 503, 504].includes(response.status))
      ) {
        throw new Error(
          `Parity download failed for ${url}: ${error instanceof Error ? error.message : String(error)}`,
          { cause: error },
        )
      }
      await sleep(1000 * 2 ** attempt)
    }
  }
}

// The same immutable, checksum-verified assets serve local previews and deployment.
export async function fetch_parity_assets(
  site_dir = fileURLToPath(new URL(`..`, import.meta.url)),
): Promise<number> {
  const pending: { file: string; entry: Asset }[] = []
  for (const kind of [`energy`, `kappa`]) {
    const manifest: Manifest = JSON.parse(
      await readFile(`${site_dir}/src/lib/parity/${kind}-parity-manifest.json`, `utf8`),
    )
    const entries = [
      manifest.base,
      ...Object.values(manifest.model_assets).flatMap(Object.values),
      ...(manifest.structure_bundles ?? []),
    ]
    const asset_dir = `${site_dir}/static/${kind}-parity/assets`
    await mkdir(asset_dir, { recursive: true })
    for (const entry of entries) {
      if (
        !entry.asset.startsWith(`${manifest.asset_prefix}-`) ||
        !/^[A-Za-z0-9._-]+-[0-9a-f]{16}\.json\.gz$/.test(entry.asset) ||
        !/^[0-9a-f]{64}$/.test(entry.sha256)
      ) {
        throw new Error(`Invalid ${kind} parity asset: ${JSON.stringify(entry)}`)
      }
      const file = `${asset_dir}/${entry.asset}`
      const existing = await readFile(file).catch((error: unknown) => {
        if (!(error instanceof Error && `code` in error && error.code === `ENOENT`)) {
          throw error
        }
        return null
      })
      if (!existing || sha256(existing) !== entry.sha256) pending.push({ file, entry })
    }
  }
  let next_idx = 0
  await Promise.all(
    Array.from({ length: Math.min(8, pending.length) }, async () => {
      while (next_idx < pending.length) {
        const { file, entry } = pending[next_idx++]
        const url = `https://github.com/janosh/matbench-discovery/releases/download/v1.0.0/${entry.asset}`
        const bytes = await download(url)
        const actual_hash = sha256(bytes)
        if (actual_hash !== entry.sha256)
          throw new Error(
            `Checksum mismatch for ${url}: expected ${entry.sha256}, received ${actual_hash}`,
          )
        const temporary_file = `${file}.${process.pid}.tmp`
        await writeFile(temporary_file, bytes)
        await rename(temporary_file, file)
      }
    }),
  )
  return pending.length
}

if (import.meta.main) {
  console.info(
    `Downloaded ${await fetch_parity_assets()} parity assets; all manifest checksums verified.`,
  )
}
