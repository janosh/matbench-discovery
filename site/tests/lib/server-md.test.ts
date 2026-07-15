import { read_md_per_system } from '$lib/server/md'
import { mkdtemp, writeFile } from 'node:fs/promises'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { gzipSync } from 'node:zlib'
import { describe, expect, it } from 'vitest'

// mirrors the pred_file CSVs written by scripts/evals/md.py: one row per DynaMat
// system, RMSEs in eV, empty cells for unavailable values (e.g. stress-less pressure)
const csv_fixture = `system,temperature_kelvin,vdos_error,energy_rmse,force_rmse,pressure_error,hardware
CsSnI3_500K_Ivor_VASP,500,3.409,0.0021,0.04175,,NVIDIA H200
bulkCu_1000K,1000,12.3,0.0005,0.02,45.6,NVIDIA H200`

describe(`read_md_per_system`, () => {
  it.each([
    [`plain`, false],
    [`gzipped`, true],
  ])(`parses a %s per-system CSV into typed rows (RMSEs in meV)`, async (_name, gzip) => {
    const tmp_dir = await mkdtemp(join(tmpdir(), `md-per-system-`))
    const csv_path = join(tmp_dir, gzip ? `metrics.csv.gz` : `metrics.csv`)
    await writeFile(csv_path, gzip ? gzipSync(csv_fixture) : csv_fixture)

    const rows = await read_md_per_system({
      metrics: { md: { pred_file: { name: csv_path } } },
    })
    if (!rows) throw new Error(`expected rows from fixture CSV`)
    expect(rows).toHaveLength(2)
    const first = rows.find((row) => row.system === `CsSnI3_500K_Ivor_VASP`)
    expect(first?.temperature_kelvin).toBe(500)
    expect(first?.vdos_error).toBeCloseTo(3.409, 3)
    // stored as 0.04175 eV/Å in the CSV, reported in meV to match the MD table
    expect(first?.force_rmse).toBeCloseTo(41.75, 2)
    expect(first?.energy_rmse).toBeCloseTo(2.1, 2)
    // empty cell -> key absent (not 0 or NaN); non-numeric cell stays a string
    expect(first).not.toHaveProperty(`pressure_error`)
    expect(first?.hardware).toBe(`NVIDIA H200`)
  })

  it.each([
    [`no pred_file`, {}],
    [`missing file`, { metrics: { md: { pred_file: { name: `missing.csv.gz` } } } }],
  ])(`returns null for %s`, async (_name, model) => {
    expect(await read_md_per_system(model)).toBeNull()
  })
})
