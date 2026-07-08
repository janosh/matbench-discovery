import { read_md_per_system } from '$lib/server/md'
import type { ModelData } from '$lib/types'
import { describe, expect, it } from 'vitest'

const model_with = (md: object | string | null): Pick<ModelData, `metrics`> =>
  ({ metrics: { md } }) as Pick<ModelData, `metrics`>

describe(`read_md_per_system`, () => {
  it(`parses a committed per-system CSV into typed rows (RMSEs in meV)`, async () => {
    const rows = await read_md_per_system(
      model_with({
        pred_file: `models/grace/grace-2l-oam/2026-07-04-grace_2l_oam-md-metrics.csv.gz`,
      }),
    )
    if (!rows) throw new Error(`expected rows from committed grace CSV`)
    expect(rows).toHaveLength(17)
    const first = rows.find((row) => row.system === `CsSnI3_500K_Ivor_VASP`)
    expect(first?.temperature_kelvin).toBe(500)
    expect(first?.vdos_error).toBeCloseTo(3.409, 3)
    // stored as 0.0417 eV/Å in the CSV, reported in meV to match the MD table
    expect(first?.force_rmse).toBeCloseTo(41.75, 2)
  })

  it.each([
    [`no md metrics`, null],
    [`md not an object`, `not available`],
    [`no pred_file`, {}],
    [`missing file`, { pred_file: `models/does/not/exist.csv.gz` }],
  ])(`returns null for %s`, async (_name, md) => {
    expect(await read_md_per_system(model_with(md))).toBeNull()
  })
})
