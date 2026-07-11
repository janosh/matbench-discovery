import { DATASETS, MODELS } from '$lib'
import { ALL_METRICS, DIATOMICS_METRICS, METADATA_COLS } from '$lib/labels'
import {
  assemble_row_data,
  format_train_set,
  metric_better_as,
  sort_models,
  targets_tooltips,
} from '$lib/metrics'
import type { TargetType } from '$lib/schema/model'
import type { ModelData } from '$lib/types'
import { describe, expect, it } from 'vitest'

describe(`targets_tooltips`, () => {
  it.each([
    [`E`, `Energy`],
    [`EF_G`, `Energy with gradient-based forces`],
    [`EFS_DM`, `Energy with direct forces, stress, and magmoms`],
    [`EF_D`, `Energy with direct forces`],
    [`EFS_G`, `Energy with gradient-based forces and stress`],
    [`EFSH_G`, `Energy with gradient-based forces, stress, and Hessian`],
  ])(`contains tooltip for %s target type`, (target, expected) => {
    expect(targets_tooltips[target as TargetType]).toBe(expected)
  })

  it(`contains all expected tooltip keys`, () => {
    expect(Object.keys(targets_tooltips)).toHaveLength(8)
  })
})

describe(`metric_better_as`, () => {
  // guards metric orientation in modeling-tasks.yml, e.g. CMDS/combined_score
  // being a score (higher=better) and errors like MAE/RMSE being lower=better
  it.each([
    [`combined_score`, `higher`],
    [`CMDS`, `higher`],
    [`vdos_error`, `lower`],
    [`F1`, `higher`],
    [`Precision`, `higher`],
    [`MAE`, `lower`],
    [`RMSE`, `lower`],
    [`nonexistent_metric`, null],
  ])(`maps %s -> %s`, (metric, expected) => {
    expect(metric_better_as(metric)).toBe(expected)
  })
})

describe(`format_train_set`, () => {
  // Get actual keys from DATASETS to use in tests
  const dataset_keys = Object.keys(DATASETS)
  const mp2022_key = dataset_keys.find((key) => key.includes(`MP 2022`))
  if (!mp2022_key) throw new Error(`No MP 2022 key found in DATASETS`)
  const mptrj_key = dataset_keys.find((key) => key.includes(`MPtrj`))
  if (!mptrj_key) throw new Error(`No MPtrj key found in DATASETS`)

  const mp2022 = DATASETS[mp2022_key]

  it(`formats single training set correctly`, () => {
    const mock_model = {
      n_training_structures: mp2022.n_structures,
      n_training_materials: mp2022.n_materials,
    }
    const result = format_train_set([mp2022_key], mock_model as ModelData)

    // Check that the result contains key information without hardcoding values
    expect(result).toContain(
      `data-sort-value="${mp2022.n_materials ?? mp2022.n_structures}"`,
    )
    expect(result).toContain(mp2022_key)
    expect(result).toContain(`materials in training set`)
  })

  it(`renders _x dataset-key suffixes as subscripts`, () => {
    const result = format_train_set([`MDR-MP PBE ω_q`], {} as ModelData)
    expect(result).toContain(`ω<sub>q</sub>`)
    expect(result).not.toContain(`ω_q`)
  })

  it(`formats multiple training sets correctly`, () => {
    const mptrj = DATASETS[mptrj_key]
    const mock_model = {
      n_training_structures: (mp2022.n_structures ?? 0) + (mptrj.n_structures ?? 0),
      n_training_materials: (mp2022.n_materials ?? 0) + (mptrj.n_materials ?? 0),
    }
    const result = format_train_set([mp2022_key, mptrj_key], mock_model as ModelData)

    // data-sort-value uses model.n_training_materials (the combined total)
    expect(result).toContain(`data-sort-value="${mock_model.n_training_materials}"`)
    expect(result).toContain(mp2022.name)
    expect(result).toContain(mptrj.name ?? mptrj_key)
  })

  it(`shows materials and structures when they differ`, () => {
    // Find a dataset with both n_materials and n_structures
    const dataset_with_both = Object.entries(DATASETS).find(
      ([_, dataset]) =>
        dataset.n_materials &&
        dataset.n_structures &&
        dataset.n_materials !== dataset.n_structures,
    )

    if (!dataset_with_both) {
      throw new Error(`No dataset with different n_materials and n_structures found`)
    }

    const [key, _dataset] = dataset_with_both
    const mock_model = {
      n_training_structures: _dataset.n_structures,
      n_training_materials: _dataset.n_materials,
    }
    const result = format_train_set([key], mock_model as ModelData)

    // Check that the result shows both materials and structures
    expect(result).toContain(`<small>(`)
    expect(result).toContain(`materials in training set (`)
    expect(result).toContain(`structures`)
  })

  it(`throws for unknown training sets`, () => {
    const mock_model = {
      n_training_structures: mp2022.n_structures,
      n_training_materials: mp2022.n_materials,
    }

    expect(() =>
      format_train_set([mp2022_key, `NonExistent`], mock_model as ModelData),
    ).toThrow(`Training set NonExistent not found in DATASETS`)
  })

  it(`formats training sets without n_materials correctly using n_structures`, () => {
    // Create a copy of a dataset and remove its n_materials property
    const mptrj = { ...DATASETS[mptrj_key] }
    const n_structures = mptrj.n_structures
    const { slug } = mptrj
    delete mptrj.n_materials

    const mock_model_struct_only = {
      n_training_structures: n_structures,
      // No n_training_materials explicitly set
    }

    // Replace the DATASETS object with our modified version
    Object.defineProperty(DATASETS, `Modified_MPtrj`, {
      value: mptrj,
      configurable: true,
    })

    try {
      const result = format_train_set(
        [`Modified_MPtrj`],
        mock_model_struct_only as ModelData,
      )

      // Should use n_structures as data-sort-value
      expect(result).toContain(`data-sort-value="${n_structures}"`)

      // Should include internal link to dataset slug page
      expect(result).toContain(`<a href="/data/${slug}"`)
    } finally {
      // Clean up by restoring original DATASETS
      delete DATASETS.Modified_MPtrj
    }
  })
})

describe(`assemble_row_data`, () => {
  const test_model_keys = [`mace-mp-0`, `chgnet-0.3.0`]
  const model_filter = (model: ModelData): boolean =>
    test_model_keys.includes(model.model_key ?? ``)
  const get_test_rows = () => assemble_row_data(`unique_prototypes`, model_filter)
  const tece_model = MODELS.find((model) => model.model_key === `tece-oam-rra-1.0`)
  if (!tece_model) throw new Error(`missing TECE-OAM-RRA-1.0 test fixture`)

  it(`returns formatted rows for selected models with expected properties`, () => {
    const rows = get_test_rows()

    expect(rows).toHaveLength(test_model_keys.length)
    const mace_row = rows.find((row) => row.Model.includes(`mace-mp-0`))
    const chgnet_row = rows.find((row) => row.Model.includes(`chgnet-0.3.0`))

    expect(mace_row?.Model).toContain(`mace-mp-0`)
    expect(mace_row?.graph_construction_radius).toBe(
      `<span data-sort-value="6">6 Å</span>`,
    )
    // n_layers should be present as either a sortable span or 'n/a'
    const n_layers_val = mace_row?.n_layers as string
    expect(n_layers_val).toMatch(/^(?:<span data-sort-value="\d+">\d+<\/span>|n\/a)$/)
    expect(chgnet_row?.Model).toContain(`chgnet-0.3.0`)
    const cps_vals = rows.map((row) => row.CPS) as number[]
    expect(cps_vals).toStrictEqual(
      cps_vals.toSorted((score_1, score_2) => score_2 - score_1),
    )
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
    { task: `diatomics`, multiplier_key: `diatomics_time_multiplier` },
    { task: `md`, multiplier_key: `md_time_multiplier` },
  ] as const)(
    `computes $task runtime multipliers relative to fastest shown model`,
    ({ task, multiplier_key }) => {
      const rows = assemble_row_data(
        `unique_prototypes`,
        (model) => model.model_key?.startsWith(`${task}-time-`) ?? false,
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

describe(`METADATA_COLS`, () => {
  it(`contains all expected metadata columns in the right order`, () => {
    const expected_labels = [
      `Model`,
      `Training Set`,
      `Targets`,
      `Date Added`,
      `Links`,
      `Training Materials`,
      `Training Structures`,
      `Ckpt License`,
      `Code License`,
      `Run Time`,
      `Org`,
    ]

    expect(Object.values(METADATA_COLS).map((col) => col.label)).toStrictEqual(
      expected_labels,
    )
  })

  it(`has the correct properties for each column`, () => {
    // Test special properties of columns
    const model_col = METADATA_COLS.model_name
    expect(model_col?.sticky).toBe(true)
    expect(model_col?.sortable).toBe(true)

    const links_col = METADATA_COLS.links
    expect(links_col?.sortable).toBe(false)
  })
})

describe(`Model Sorting Logic`, () => {
  // Create a set of test models that can be reused across multiple test cases
  const create_test_models = () =>
    [
      {
        model_name: `AAA Model`,
        model_key: `aaa_model`,
        metrics: {
          discovery: {
            unique_prototypes: { F1: 0.9, Accuracy: 0.85, missing_preds: 0 },
            pred_col: `is_stable`,
          },
          phonons: { kappa_103: { κ_SRME: 0.9 } },
        },
      },
      {
        model_name: `MMM Model`,
        model_key: `mmm_model`,
        metrics: {
          discovery: {
            unique_prototypes: { F1: 0.7, Accuracy: NaN, missing_preds: 2 }, // NaN Accuracy + non-zero missing_preds
            pred_col: `is_stable`,
          },
          phonons: { kappa_103: { κ_SRME: 0.5 } },
        },
      },
      {
        model_name: `ZZZ Model`,
        model_key: `zzz_model`,
        metrics: {
          discovery: {
            unique_prototypes: { F1: 0.5, Accuracy: 0.6, missing_preds: 5 },
            pred_col: `is_stable`,
          },
          phonons: { kappa_103: { κ_SRME: 0.2 } },
        },
      },
      {
        model_name: `Missing Data Model`,
        model_key: `missing_model`,
        metrics: {
          discovery: { unique_prototypes: { F1: 0.4 } },
          phonons: `not applicable` as const,
        },
      },
    ] as unknown as ModelData[]

  // Create test models with edge cases
  const create_edge_case_models = () =>
    [
      // Model with completely missing metrics
      {
        model_name: `No Metrics Model`,
        model_key: `no_metrics_model`,
      },
      // Model with empty metrics object
      {
        model_name: `Empty Metrics Model`,
        model_key: `empty_metrics_model`,
        metrics: {},
      },
      // Model with extreme values
      {
        model_name: `Extreme Values Model`,
        model_key: `extreme_model`,
        metrics: {
          discovery: {
            unique_prototypes: {
              F1: Number.MAX_VALUE,
              Accuracy: Number.MIN_VALUE,
              R2: -Infinity,
              RMSE: Infinity,
            },
          },
          phonons: {
            kappa_103: { κ_SRME: 0 },
          },
        },
      },
      // Model with all undefined values for metrics
      {
        model_name: `Undefined Metrics Model`,
        model_key: `undefined_model`,
        metrics: {
          discovery: {
            unique_prototypes: {
              F1: undefined,
              Accuracy: undefined,
            },
          },
          phonons: {
            kappa_103: { κ_SRME: undefined },
          },
        },
      },
    ] as unknown as ModelData[]

  // mmm_model has NaN Accuracy and missing_model lacks metrics, so these cases also
  // assert NaN/missing values sort last for every metric and direction (incl. symmetric
  // Accuracy asc/desc)
  it.each<[string, `asc` | `desc`, string[]]>([
    [
      `${ALL_METRICS.F1.path}.${ALL_METRICS.F1.key}`,
      `desc`,
      [`aaa_model`, `mmm_model`, `zzz_model`, `missing_model`],
    ],
    [
      `${ALL_METRICS.Accuracy.path}.${ALL_METRICS.Accuracy.key}`,
      `asc`,
      [`zzz_model`, `aaa_model`, `mmm_model`, `missing_model`],
    ],
    [
      `${ALL_METRICS.Accuracy.path}.${ALL_METRICS.Accuracy.key}`,
      `desc`,
      [`aaa_model`, `zzz_model`, `mmm_model`, `missing_model`],
    ],
    [
      `${ALL_METRICS.κ_SRME.path}.${ALL_METRICS.κ_SRME.key}`,
      `asc`,
      [`zzz_model`, `mmm_model`, `aaa_model`, `missing_model`],
    ],
    [
      `metrics.discovery.unique_prototypes.missing_preds`,
      `asc`,
      [`aaa_model`, `mmm_model`, `zzz_model`, `missing_model`],
    ],
  ])(`sorts test models by %s (%s)`, (metric, order, expected_order) => {
    const sorted = create_test_models().toSorted(sort_models(metric, order))
    expect(sorted.map((model) => model.model_key)).toStrictEqual(expected_order)
  })

  it(`returns 0 for two models that both have NaN values (consistent comparator)`, () => {
    // regression: NaN-vs-NaN returned 1 both ways, breaking antisymmetry -> order-dependent sort
    const { Accuracy } = ALL_METRICS
    const sort_by = `${Accuracy.path}.${Accuracy.key}`
    const make_nan_model = (model_key: string) =>
      ({
        model_key,
        metrics: { discovery: { unique_prototypes: { Accuracy: NaN } } },
      }) as unknown as ModelData
    const [nan_1, nan_2] = [make_nan_model(`nan_1`), make_nan_model(`nan_2`)]
    const valid = {
      model_key: `valid`,
      metrics: { discovery: { unique_prototypes: { Accuracy: 0.5 } } },
    } as unknown as ModelData

    const compare = sort_models(sort_by, `desc`)
    expect(compare(nan_1, nan_2)).toBe(0)
    expect(compare(nan_2, nan_1)).toBe(0)
    // NaN still sorts after real values in both argument orders
    expect(compare(valid, nan_1)).toBeLessThan(0)
    expect(compare(nan_1, valid)).toBeGreaterThan(0)

    // sorting NaN-heavy arrays is now stable regardless of initial order
    const models = [nan_1, valid, nan_2]
    const fwd = models.toSorted(compare).map((model) => model.model_key)
    const rev = models
      .toReversed()
      .toSorted(compare)
      .map((model) => model.model_key)
    expect(fwd[0]).toBe(`valid`)
    expect(rev[0]).toBe(`valid`)
  })

  it(`sorts models by model_name correctly`, () => {
    const asc_keys = create_test_models()
      .toSorted(sort_models(`model_name`, `asc`))
      .map((model) => model.model_key)
    expect(asc_keys).toStrictEqual([
      `aaa_model`,
      `missing_model`,
      `mmm_model`,
      `zzz_model`,
    ])

    // descending is the exact reverse of ascending
    const desc_keys = create_test_models()
      .toSorted(sort_models(`model_name`, `desc`))
      .map((model) => model.model_key)
    expect(desc_keys).toStrictEqual(asc_keys.toReversed())
  })

  it(`handles edge cases with missing or extreme metric values`, () => {
    const edge_case_models = create_edge_case_models()
    const regular_models = create_test_models()
    const combined_models = [...regular_models, ...edge_case_models]

    // Test sorting with κ_SRME where one model has zero value
    const { κ_SRME } = ALL_METRICS
    const sort_by_path = `${κ_SRME.path}.${κ_SRME.key}`
    const sorted_by_kappa = combined_models.toSorted(sort_models(sort_by_path, `asc`))

    // Zero value should be first for asc
    expect(sorted_by_kappa[0].model_key).toBe(`extreme_model`)
  })

  it(`sorts models by runtime correctly, treating 0 as infinity`, () => {
    const models = Object.entries({ a: 10, b: 0, c: 5, d: 0 }).map(
      ([model_key, run_time]) => ({
        model_key: `model_${model_key}`,
        'Run Time': run_time,
      }),
    ) as unknown as ModelData[]

    // Test ascending sort (runtime 0 should be last)
    const sorted_asc = models.toSorted(sort_models(`Run Time`, `asc`))
    // Check the non-zero values are sorted correctly first
    expect(sorted_asc.slice(0, 2).map((model) => model.model_key)).toStrictEqual([
      `model_c`,
      `model_a`,
    ])
    // Check the zero values are at the end (order between them is not guaranteed)
    expect(
      sorted_asc
        .slice(2)
        .map((model) => model.model_key)
        .toSorted((key_1, key_2) => (key_1 ?? ``).localeCompare(key_2 ?? ``)),
    ).toStrictEqual([`model_b`, `model_d`])

    // Test descending sort (runtime 0 should be first)
    const sorted_desc = models.toSorted(sort_models(`Run Time`, `desc`))
    // Check the zero values are at the beginning (order between them is not guaranteed)
    expect(
      sorted_desc
        .slice(0, 2)
        .map((model) => model.model_key)
        .toSorted((key_1, key_2) => (key_1 ?? ``).localeCompare(key_2 ?? ``)),
    ).toStrictEqual([`model_b`, `model_d`])
    // Check the non-zero values are sorted correctly after
    expect(sorted_desc.slice(2).map((model) => model.model_key)).toStrictEqual([
      `model_a`,
      `model_c`,
    ])
  })

  it(`handles sorting with all models having missing metric`, () => {
    const all_missing_models = [
      { model_name: `Model A`, model_key: `model_a` },
      { model_name: `Model B`, model_key: `model_b` },
    ] as unknown as ModelData[]

    // Sort by a metric that none of the models have
    const sorted_models = all_missing_models.toSorted(sort_models(`F1`, `desc`))

    // Order should be preserved when all models are missing the metric
    expect(sorted_models[0].model_key).toBe(`model_a`)
    expect(sorted_models[1].model_key).toBe(`model_b`)
  })

  it(`maintains original order for equivalent values`, () => {
    const models_with_same_values = [1, 2, 3].map((idx) => ({
      model_name: `Model ${idx}`,
      model_key: `model_${idx}`,
      metrics: { discovery: { unique_prototypes: { F1: 0.8 } } },
    })) as unknown as ModelData[]

    // Original order should be preserved when sorting identical F1 values
    const sorted_models = models_with_same_values.toSorted(sort_models(`F1`, `desc`))
    const sorted_keys = sorted_models.map((model) => model.model_key)
    expect(sorted_keys).toStrictEqual([`model_1`, `model_2`, `model_3`])
  })

  it(`throws an error when sorting by unexpected type`, () => {
    const models_with_mixed_types = [
      { model_name: `Model A`, model_key: `model_a`, some_metric: 1 },
      { model_name: `Model B`, model_key: `model_b`, some_metric: `string` },
    ] as unknown as ModelData[]

    // Expect an error when trying to sort number and string
    expect(() =>
      models_with_mixed_types.sort(sort_models(`some_metric`, `desc`)),
    ).toThrow(/Unexpected type.*encountered sorting by key/)
  })

  it.each([`asc`, `desc`] as const)(`sorts nulls last in %s order`, (order) => {
    const models_with_null = [
      { model_name: `Model Null 1`, model_key: `null_1`, metric: null },
      { model_name: `Model Null 2`, model_key: `null_2`, metric: null },
      { model_name: `Model Val`, model_key: `val_1`, metric: 10 },
    ] as unknown as ModelData[]

    expect(
      models_with_null
        .toSorted(sort_models(`metric`, order))
        .map((model) => model.model_key),
    ).toStrictEqual([`val_1`, `null_1`, `null_2`])
  })
})
