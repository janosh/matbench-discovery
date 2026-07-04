import type { ModelData } from '$lib'
import {
  calculate_cmds,
  type CmdsConfig,
  DEFAULT_CMDS_CONFIG,
  update_models_cmds,
} from '$lib/md_combined_score.svelte'
import { describe, expect, it } from 'vitest'

const clone_config = (): CmdsConfig => structuredClone(DEFAULT_CMDS_CONFIG)

describe(`calculate_cmds`, () => {
  it.each([
    // equal weights reproduce the formerly stored formula 1 - mean(errors)/100
    [{ adf_error: 10, vdos_error: 20, pressure_error: 30 }, 0.8],
    [{ adf_error: 0, vdos_error: 0, pressure_error: 0 }, 1],
    [{ adf_error: 100, vdos_error: 100, pressure_error: 100 }, 0],
    // per-component subscores clamp at 0, so >100% errors can't go negative
    [{ adf_error: 150, vdos_error: 0, pressure_error: 0 }, 2 / 3],
  ])(`equal-weight score for %o = %f`, (errors, expected) => {
    expect(calculate_cmds(errors, clone_config())).toBeCloseTo(expected, 10)
  })

  it.each([
    [{ adf_error: 10, vdos_error: 20 }], // missing component
    [{ adf_error: 10, vdos_error: 20, pressure_error: Number.NaN }], // NaN component
  ])(`returns null when a non-zero-weight component is invalid: %o`, (errors) => {
    expect(calculate_cmds(errors, clone_config())).toBeNull()
  })

  it(`skips missing components with zero weight`, () => {
    const config = clone_config()
    config.pressure_error.weight = 0
    // mean of (1 - 10/100) and (1 - 20/100) weighted 1/3 each = 0.85
    expect(calculate_cmds({ adf_error: 10, vdos_error: 20 }, config)).toBeCloseTo(0.85)
  })

  it(`returns null when all weights are zero (score undefined, not worst)`, () => {
    const config = clone_config()
    for (const key of Object.keys(config) as (keyof CmdsConfig)[]) {
      config[key].weight = 0
    }
    expect(calculate_cmds({}, config)).toBeNull()
  })
})

describe(`update_models_cmds`, () => {
  it(`writes CMDS into metrics.md and deletes it when components are missing`, () => {
    const models = [
      { metrics: { md: { adf_error: 10, vdos_error: 20, pressure_error: 30 } } },
      { metrics: { md: { adf_error: 10, combined_score: 0.99 } } }, // stale score
      { metrics: {} },
    ] as unknown as ModelData[]

    update_models_cmds(models, clone_config())

    const scored_md = models[0].metrics?.md as { combined_score?: number }
    expect(scored_md.combined_score).toBeCloseTo(0.8, 10)
    // incomplete components: stale score removed, not overwritten with NaN (a NaN
    // field would count as present in scatter model counts and export as 'NaN')
    expect(models[1].metrics?.md).not.toHaveProperty(`combined_score`)
    expect(models[2].metrics?.md).toBeUndefined()
  })
})
