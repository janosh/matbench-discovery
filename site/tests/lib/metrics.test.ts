import DATASETS from '$data/datasets.yml'
import { MODELS } from '#lib/models.svelte.js'
import {
  ALL_METRICS,
  DIATOMICS_METRICS,
  MD_METRICS,
  PHONON_METRICS,
} from '#lib/labels.js'
import {
  assemble_row_data,
  format_train_set,
  missing_metric_reason,
  sort_models,
} from '#lib/metrics.js'
import type { ModelData } from '#lib/types.js'
import { describe, expect, it } from 'vite-plus/test'

// labels own their direction (no fallback list), so a metric added without `better` would
// silently lose its heatmap orientation and better-hint
it(`every metric declares a direction except the deliberate exceptions`, () => {
  const no_direction = Object.values(ALL_METRICS)
    .filter(({ better }) => better !== `higher` && better !== `lower`)
    .map(({ key }) => key)
  expect(no_direction).toEqual([
    `κ_SRD`,
    `symmetry_increase_1e-2`,
    `symmetry_increase_1e-5`,
  ])
})

describe(`format_train_set`, () => {
  const mp2022 = DATASETS[`MP 2022`]
  const mptrj = DATASETS.MPtrj

  it.each([`MP 2022`, `OC22`, `MAD-1.6`, `OCx24`])(
    `formats single training set %s correctly`,
    (key) => {
      const dataset = DATASETS[key]
      const result = format_train_set([key], {
        n_training_structures: dataset.n_structures,
        n_training_materials: dataset.n_materials,
      })

      expect(result).toContain(`/data/${dataset.slug}`)
      if (dataset.n_structures === null || dataset.n_materials == null) {
        expect(result).toContain(`size not fully reported`)
        expect(result).toContain(
          dataset.n_structures === null
            ? `unknown number of structures`
            : `material count not reported`,
        )
        expect(result).toContain(`n/a`)
        expect(result).not.toContain(`data-sort-value=`)
        return
      }
      expect(result).toContain(`data-sort-value="${dataset.n_materials}"`)
      expect(result).toContain(key)
      expect(result).toContain(`materials in training set`)
    },
  )

  it(`renders _x dataset-key suffixes as subscripts`, () => {
    const result = format_train_set([`MDR-MP PBE ω_q`], {})
    expect(result).toContain(`ω<sub>q</sub>`)
    expect(result).not.toContain(`ω_q`)
  })

  it(`formats multiple training sets correctly`, () => {
    const model = {
      n_training_structures: (mp2022.n_structures ?? 0) + (mptrj.n_structures ?? 0),
      n_training_materials: (mp2022.n_materials ?? 0) + (mptrj.n_materials ?? 0),
    }
    const result = format_train_set([`MP 2022`, `MPtrj`], model)

    expect(result).toContain(`data-sort-value="${model.n_training_materials}"`)
    expect(result).toContain(mp2022.name)
    expect(result).toContain(mptrj.name)
  })

  it(`shows materials and structures when they differ`, () => {
    const result = format_train_set([`MPtrj`], {
      n_training_structures: mptrj.n_structures,
      n_training_materials: mptrj.n_materials,
    })
    expect(result).toContain(`<small>(`)
    expect(result).toContain(`materials in training set (`)
    expect(result).toContain(`structures`)
  })

  it(`throws for unknown training sets`, () => {
    expect(() => format_train_set([`MP 2022`, `NonExistent`], {})).toThrow(
      `Training set NonExistent not found in DATASETS`,
    )
  })

  it(`never presents structures as materials when the material count is absent`, () => {
    const { n_structures, n_materials, slug } = DATASETS[`MAD-1.6`]
    expect(n_materials).toBeUndefined()
    const result = format_train_set([`MAD-1.6`], { n_training_structures: n_structures })
    expect(result).toContain(`362,646 structures; material count not reported`)
    expect(result).not.toContain(`362,646 materials`)
    expect(result).not.toContain(`0 materials`)
    expect(result).not.toContain(`data-sort-value=`)
    expect(result).toContain(`<a href="/data/${slug}"`)
  })
})

describe(`assemble_row_data`, () => {
  const test_model_keys = [`mace-mp-0`, `chgnet-0.3.0`]
  const model_filter = (model: ModelData): boolean =>
    test_model_keys.includes(model.model_key)
  const get_test_rows = () => assemble_row_data(`unique_prototypes`, model_filter)
  const tece_model = MODELS.find((model) => model.model_key === `tece-oam-rra-1.0`)
  if (!tece_model) throw new Error(`missing TECE-OAM-RRA-1.0 test fixture`)

  it(`returns formatted rows for selected models with expected properties`, () => {
    const rows = get_test_rows()

    expect(rows).toHaveLength(test_model_keys.length)
    const mace_row = rows.find((row) => row.Model.includes(`mace-mp-0`))
    const chgnet_row = rows.find((row) => row.Model.includes(`chgnet-0.3.0`))

    expect(mace_row?.model_key).toBe(`mace-mp-0`)
    expect(mace_row?.graph_construction_radius).toBe(
      `<span data-sort-value="6">6 Å</span>`,
    )
    // Missing metadata carries its own explanation.
    const n_layers_val = mace_row?.n_layers as string
    expect(n_layers_val).toMatch(
      /^<span (?:data-sort-value="\d+">\d+|data-title="[^"]+">n\/a)<\/span>$/,
    )
    expect(chgnet_row?.Model).toContain(`chgnet-0.3.0`)
    for (const [metric, expected] of [
      [PHONON_METRICS.κ_SRME, 0.6823],
      [PHONON_METRICS.κ_SRE, 0.471],
      [PHONON_METRICS.κ_SRD, -0.2845],
      [PHONON_METRICS.κ_failure_rate, 0.0291],
      [PHONON_METRICS.imaginary_mode_rate, 0],
      [PHONON_METRICS.spectrum_w1, 0.8709],
    ] as const) {
      expect(mace_row?.[metric.key]).toBe(expected)
    }
    const cps_vals = rows.map((row) => row.CPS) as number[]
    expect(cps_vals).toStrictEqual(
      cps_vals.toSorted((score_1, score_2) => score_2 - score_1),
    )
    expect(JSON.stringify(mace_row?.Links)).not.toMatch(/icon|title|viewBox/)
  })

  it(`includes task-only models without discovery metrics`, () => {
    const model_key = `task-only-regression`
    const task_only_model = {
      ...tece_model,
      model_key,
      model_name: `Task-only regression`,
      metrics: { diatomics: { pbe_energy_mae: 1 } },
    } as ModelData

    const rows = assemble_row_data(
      `unique_prototypes`,
      (model) => model.model_key === model_key,
      () => true,
      [task_only_model],
    )
    expect(rows).toHaveLength(1)
    expect(rows[0].F1).toBeUndefined()
    expect(rows[0][DIATOMICS_METRICS.pbe_energy_mae.key]).toBe(1)
  })

  it.each([
    [ALL_METRICS.RMSD, {}, `predicts only energies`, false, `E`],
    [PHONON_METRICS.κ_SRME, {}, `requires forces`, false, `E`],
    [MD_METRICS.md_combined_score, {}, `no results reported`, true],
    [DIATOMICS_METRICS.diatomics_combined_score, {}, `no results reported`, true],
    [
      ALL_METRICS.RMSD,
      {
        geo_opt: { status: `not_applicable`, reason: `No structure relaxation support.` },
      },
      `unsupported. No structure relaxation support.`,
      false,
    ],
    [
      PHONON_METRICS.κ_SRME,
      { phonons: { status: `pending`, reason: `Evaluation is queued.` } },
      `pending. Evaluation is queued.`,
      true,
    ],
    [
      PHONON_METRICS.κ_SRME,
      { phonons: { status: `not_available`, reason: `Predictions were not retained.` } },
      `unavailable. Predictions were not retained.`,
      true,
    ],
    [
      PHONON_METRICS.κ_SRME,
      {
        phonons: {
          status: `partial`,
          reason: `Missing conductivity predictions.`,
          kappa_103: {},
        },
      },
      `incomplete results. Missing conductivity predictions.`,
      true,
    ],
    [
      MD_METRICS.md_combined_score,
      { md: { adf_error: 0, vdos_error: 0, pressure_error: 0 } },
      `runtime not reported`,
      false,
    ],
    [
      MD_METRICS.md_max_gpu_mem_gb,
      { md: { vdos_error: 0 } },
      `peak memory not reported`,
      false,
    ],
    [ALL_METRICS.CPS, {}, `requires forces`, true, `E`],
    [ALL_METRICS.CPS, {}, `no results reported`, true],
  ] as const)(
    `explains missing results (%#)`,
    (label, metrics, expected, invite, targets: ModelData['targets'] = `EFS_G`) => {
      const model = { ...tece_model, targets, metrics } as ModelData
      const reason = missing_metric_reason(model, label)
      expect(reason).toContain(expected)
      expect(
        reason.match(/Contributions welcome to add missing model predictions/g) ?? [],
      ).toHaveLength(invite ? 1 : 0)
      if (invite)
        expect(reason).toMatch(
          /Contributions welcome to add missing model predictions\.$/,
        )
    },
  )

  it.each([
    ALL_METRICS.CPS,
    ALL_METRICS.RMSD,
    PHONON_METRICS.κ_SRME,
    MD_METRICS.md_combined_score,
    MD_METRICS.md_run_time_sec,
    MD_METRICS.md_max_gpu_mem_gb,
    DIATOMICS_METRICS.diatomics_combined_score,
  ])(`explains GNoME's missing $label without soliciting predictions`, (label) => {
    const gnome = MODELS.find(({ model_key }) => model_key === `gnome`)
    if (!gnome) throw new Error(`Missing GNoME test fixture`)
    const reason = missing_metric_reason(gnome, label)
    expect(reason).not.toMatch(/Contributions welcome|not evaluated yet/)
    expect(reason).toMatch(/Model weights are not publicly available\.$/)
    expect(reason.match(/Model weights are not publicly available\./g)).toHaveLength(1)
  })

  it.each([
    { task: `diatomics`, multiplier_key: `diatomics_time_multiplier` },
    { task: `md`, multiplier_key: `md_time_multiplier` },
  ] as const)(
    `computes $task runtime multipliers relative to fastest shown model`,
    ({ task, multiplier_key }) => {
      const rows = assemble_row_data(
        `unique_prototypes`,
        (model) => model.model_key.startsWith(`${task}-time-`),
        () => true,
        Object.entries({
          Fast: 10,
          Medium: 20,
          Slow: 40,
          Zero: 0,
          Missing: undefined,
          Infinite: Infinity,
          NaN: Number.NaN,
        }).map(([model_name, run_time_sec]) => ({
          ...tece_model,
          model_key: `${task}-time-${model_name.toLowerCase()}`,
          model_name,
          metrics: {
            ...tece_model.metrics,
            [task]: run_time_sec === undefined ? {} : { run_time_sec },
          },
        })),
      )

      expect(
        Object.fromEntries(
          rows.map((row) => [
            row.model_name,
            (row as Record<string, unknown>)[multiplier_key],
          ]),
        ),
      ).toEqual({
        Fast: 1,
        Medium: 2,
        Slow: 4,
        Zero: undefined,
        Missing: undefined,
        Infinite: undefined,
        NaN: undefined,
      })
    },
  )

  it.each([
    {
      model_key: `sevennet-l3i5`,
      diatomics: undefined,
      expected_title: `Diatomics metrics exclude He-He due to exploding errors`,
    },
    {
      model_key: `mixed-diatomics-exclusions`,
      diatomics: {
        excluded_formula_reasons: {
          'H-H': `unsupported "quoted" reason`,
          'He-He': `exploding errors`,
          'Li-Li': `exploding errors`,
        },
      },
      expected_title:
        `Diatomics metrics exclude H-H due to unsupported &quot;quoted&quot; reason; ` +
        `He-He, Li-Li due to exploding errors`,
    },
  ])(
    `renders reason-aware diatomics exclusion tooltip for $model_key`,
    ({ model_key, diatomics, expected_title }) => {
      const test_models = [
        ...MODELS,
        ...(diatomics
          ? [
              {
                ...tece_model,
                model_key,
                model_name: `Mixed Diatomics Exclusions`,
                metrics: { ...tece_model.metrics, diatomics },
              } as ModelData,
            ]
          : []),
      ]
      const [row] = assemble_row_data(
        `unique_prototypes`,
        (model) => model.model_key === model_key,
        () => true,
        test_models,
      )

      expect(row?.Model).toContain(`title="${expected_title}"`)
      expect(row?.Model).toContain(`aria-label="${expected_title}"`)
    },
  )
})

// ordering rules are matterviz's sort_table_rows; this covers the site-side mapping
it.each([
  [`Model`, `asc`, [`a9`, `a10`, `b`, `c`]],
  [`Model`, `desc`, [`c`, `b`, `a10`, `a9`]],
  [`metrics.F1`, `asc`, [`b`, `a10`, `a9`, `c`]],
  [`metrics.F1`, `desc`, [`a10`, `b`, `a9`, `c`]],
] as const)(`sort_models by %s %s`, (sort_by, order, expected) => {
  const models = [
    { model_name: `b`, metrics: { F1: 0.2 } },
    { model_name: `a10`, metrics: { F1: 0.9 } },
    { model_name: `a9`, metrics: { F1: NaN } },
    { model_name: `c`, metrics: {} },
  ] as unknown as ModelData[]
  const names = sort_models(models, sort_by, order).map(({ model_name }) => model_name)
  expect(names).toStrictEqual(expected)
})
