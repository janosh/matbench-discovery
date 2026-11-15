import kappa_analysis from '#figs/kappa-103-analysis.jsonl'
import {
  kappa_material,
  kappa_material_summaries,
  spacegroup_html,
} from '#lib/materials.js'
import { describe, expect, it } from 'vite-plus/test'

describe(`kappa_material`, () => {
  it(`slices one material's per-model rows, sorted by SRME with missing models last`, () => {
    const [first_model, second_model, ...rest] = kappa_analysis.models
    if (!first_model || !second_model) throw new Error(`need two models`)
    const idx = 0
    const material_id = kappa_analysis.material_ids[idx]
    const kappa_dft = kappa_analysis.kappa_dft[idx]
    if (kappa_dft === null) throw new Error(`fixture material has no DFT κ`)
    const with_values = (model: typeof first_model, values: Record<string, unknown>) => ({
      ...model,
      ...Object.fromEntries(
        Object.entries(values).map(([key, value]) => {
          const column = [...(model[key as keyof typeof model] as unknown[])]
          column[idx] = value
          return [key, column]
        }),
      ),
    })
    // first model: censored failure; second model: no prediction at all
    const analysis = {
      ...kappa_analysis,
      models: [
        with_values(first_model, {
          srme: 2,
          srme_censored: true,
          imag_modes: true,
          kappa_ml: kappa_dft / 4,
        }),
        with_values(second_model, { srme: null, kappa_ml: null, srme_censored: null }),
        ...rest,
      ],
    }
    const material = kappa_material(analysis, material_id)

    expect(material.rows).toHaveLength(kappa_analysis.models.length)
    const srmes = material.rows.map((row) => row.srme)
    expect(srmes.at(-1)).toBeNull() // unscored model sorts last
    const scored = srmes.filter((srme) => srme !== null)
    expect(scored).toEqual(scored.toSorted((srme_a, srme_b) => srme_a - srme_b))

    const failed = material.rows.find((row) => row.model_key === first_model.model_key)
    expect(failed).toMatchObject({ srme: 2, srme_censored: true, imag_modes: true })
    expect(failed?.kappa_ratio).toBeCloseTo(0.25, 12)
    expect(material.n_models).toBe(kappa_analysis.models.length - 1)
    expect(material.worst?.srme).toBe(2)
    expect(material.best?.srme).toBe(scored[0])
    expect([material.prev_id, material.next_id]).toEqual([
      null,
      kappa_analysis.material_ids[1],
    ])
  })

  it(`throws for materials outside the κ-103 set`, () => {
    expect(() => kappa_material(kappa_analysis, `mp-0`)).toThrow(`mp-0 not in the κ-103`)
  })
})

describe(`kappa_material median rank`, () => {
  it(`ranks median SRME from lowest (1) to highest across all materials`, () => {
    const summaries = kappa_material_summaries(kappa_analysis)
    const medians = summaries.flatMap(({ median_srme }) => median_srme ?? [])
    const ranked = summaries.map(({ material_id }) =>
      kappa_material(kappa_analysis, material_id),
    )
    for (const { median_srme_rank, n_ranked } of ranked) {
      expect(n_ranked).toBe(medians.length)
      expect(median_srme_rank).toBeGreaterThanOrEqual(1)
      expect(median_srme_rank).toBeLessThanOrEqual(n_ranked)
    }
    const by_median = ranked.toSorted(
      (mat_a, mat_b) => (mat_a.median_srme ?? 0) - (mat_b.median_srme ?? 0),
    )
    expect(by_median[0]?.median_srme_rank).toBe(1)
    expect(by_median.at(-1)?.median_srme_rank).toBe(medians.length)
  })
})

describe(`kappa_material_summaries`, () => {
  it(`summarizes every material, with median SRME between best and worst`, () => {
    const summaries = kappa_material_summaries(kappa_analysis)
    expect(summaries.map((summary) => summary.material_id)).toEqual(
      kappa_analysis.material_ids,
    )
    for (const { best, worst, median_srme } of summaries) {
      if (!best || !worst || median_srme === null) continue
      expect(best.srme).toBeLessThanOrEqual(median_srme)
      expect(median_srme).toBeLessThanOrEqual(worst.srme)
    }
  })
})

describe(`spacegroup_html`, () => {
  it.each([
    [186, `P6<sub>3</sub>mc (186)`],
    [225, `Fm<span style="text-decoration: overline">3</span>m (225)`],
    [1, `P1 (1)`],
    [999, `999`],
  ])(`spacegroup %i -> %s`, (spg_num, expected) => {
    expect(spacegroup_html(spg_num)).toBe(expected)
  })
})
