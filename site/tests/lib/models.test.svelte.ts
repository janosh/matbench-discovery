import {
  calculate_cps,
  CPS_CONFIG,
  DEFAULT_CPS_CONFIG,
} from '$lib/combined-scores.svelte'
import { attach_style, order_models } from '$lib/fig-helpers'
import { ALL_METRICS } from '$lib/labels'
import {
  ALL_TRAINING_SETS,
  calculate_training_sizes,
  find_best_model,
  make_table_filters,
  MODEL_METADATA_PATHS,
  MODELS,
  update_models_cps,
} from '$lib/models.svelte'
import type { ModelData } from '$lib/types'
import { per_element_each_errors as per_elem_each_errors } from '$lib/per-element-errors'
import { load as yaml_load } from 'js-yaml'
import { readdirSync, readFileSync } from 'node:fs'
import path from 'node:path'
import { describe, expect, it } from 'vitest'

describe(`calculate_training_sizes`, () => {
  it(`should return zeros for empty training sets`, () => {
    const result = calculate_training_sizes([])
    expect(result.total_materials).toBe(0)
    expect(result.total_structures).toBe(0)
  })

  it(`should calculate totals correctly for a single dataset`, () => {
    const result = calculate_training_sizes([`MP 2022`])
    expect(result.total_materials).toBe(154_719)
    expect(result.total_structures).toBe(154_719)
  })

  it(`should calculate totals correctly for multiple datasets`, () => {
    const result = calculate_training_sizes([`MP 2022`, `MPtrj`])
    expect(result.total_materials).toBe(300_642)
    expect(result.total_structures).toBe(1_735_114)
  })

  it(`should use n_structures as n_materials when n_materials is not specified`, () => {
    const result = calculate_training_sizes([`MPF`])
    expect(result.total_materials).toBe(62_783) // Should equal n_structures
    expect(result.total_structures).toBe(188_349)
  })

  it(`throws for unknown training sets`, () => {
    expect(() => calculate_training_sizes([`NonExistentDataset`, `MP 2022`])).toThrow(
      `Training set NonExistentDataset not found in DATASETS`,
    )
  })
})

describe(`MODELS array`, () => {
  it(`should be a non-empty array`, () => {
    expect(Array.isArray(MODELS)).toBe(true)
    expect(MODELS.length).toBeGreaterThan(0)
  })

  it(`should have processed models with calculated properties`, () => {
    // Skip if no models available in test environment
    if (MODELS.length === 0) return

    // Check that models have the expected structure
    const model = MODELS[0]
    expect(model).toHaveProperty(`dirname`)
    expect(model).toHaveProperty(`metadata_file`)
    expect(model).toHaveProperty(`color`)
    expect(model).toHaveProperty(`n_training_materials`)
    expect(model).toHaveProperty(`n_training_structures`)
  })

  it(`includes logos from non-lead author affiliations`, () => {
    const mirror_physics_model = MODELS.find((model) =>
      model.authors?.some((author) => author.affiliation === `Mirror Physics`),
    )
    if (!mirror_physics_model) throw new Error(`missing Mirror Physics model`)

    expect(mirror_physics_model.authors[0]?.affiliation).not.toBe(`Mirror Physics`)
    expect(mirror_physics_model.org_logos?.map((logo) => logo.name)).toContain(
      `Mirror Physics`,
    )
  })
})

describe(`MODEL_METADATA_PATHS`, () => {
  it(`should be defined and be an object`, () => {
    expect(typeof MODEL_METADATA_PATHS).toBe(`object`)
    const model_keys = new Set(
      Object.values(MODEL_METADATA_PATHS).map((model) => model.model_key),
    )
    const data_keys = Object.keys(per_elem_each_errors).filter(
      (key) => ![`MP Occurrences`, `Test set standard deviation`].includes(key),
    )
    expect(data_keys.every((key) => model_keys.has(key))).toBe(true)
    expect(per_elem_each_errors).toHaveProperty(`mace-mp-0`)
    expect(per_elem_each_errors).not.toHaveProperty(`MACE-MP-0`)
  })
})

describe(`make_table_filters`, () => {
  it(`default filters only hide energy-only models`, () => {
    const filters = make_table_filters()
    expect(filters.n_active).toBe(0) // the default require-F isn't counted as active
    for (const model of MODELS) {
      expect(filters.matches(model), model.model_name).toBe(model.targets !== `E`)
    }
  })

  it.each<{
    training: Record<string, `require` | `exclude`>
    training_set: string[]
    expected: boolean
  }>([
    // require: only models trained on the dataset
    { training: { MPtrj: `require` }, training_set: [`MPtrj`], expected: true },
    { training: { MPtrj: `require` }, training_set: [`OMat24`], expected: false },
    // exclude: only models NOT trained on the dataset
    { training: { OMat24: `exclude` }, training_set: [`MPtrj`], expected: true },
    {
      training: { OMat24: `exclude` },
      training_set: [`MPtrj`, `OMat24`],
      expected: false,
    },
    // multiple requires AND together
    {
      training: { MPtrj: `require`, sAlex: `require` },
      training_set: [`MPtrj`, `sAlex`, `OMat24`],
      expected: true,
    },
    {
      training: { MPtrj: `require`, sAlex: `require` },
      training_set: [`MPtrj`],
      expected: false,
    },
  ])(
    `training filter $training matches $training_set -> $expected`,
    ({ training, training_set, expected }) => {
      const filters = make_table_filters()
      filters.training = training
      expect(filters.matches({ training_set, targets: `EFS_G` })).toBe(expected)
    },
  )

  it.each([
    { openness: [`OSOD`], model_openness: `OSOD`, expected: true },
    { openness: [`OSOD`], model_openness: `CSCD`, expected: false },
    { openness: [`OSOD`], model_openness: undefined, expected: true }, // defaults OSOD
    { openness: [`OSCD`, `CSCD`], model_openness: `OSOD`, expected: false },
  ] as const)(
    `openness filter $openness matches $model_openness -> $expected`,
    ({ openness, model_openness, expected }) => {
      const filters = make_table_filters()
      filters.openness = [...openness]
      expect(
        filters.matches({
          training_set: [`MPtrj`],
          openness: model_openness,
          targets: `EFS_G`,
        }),
      ).toBe(expected)
    },
  )

  it.each([
    // default targets filter (require F) hides energy-only models
    { model_targets: `E`, expected: false },
    { model_targets: `EFS_G`, expected: true },
    // exclude magmoms
    { targets: { F: `require`, M: `exclude` }, model_targets: `EFS_GM`, expected: false },
    { targets: { F: `require`, M: `exclude` }, model_targets: `EFS_G`, expected: true },
    // require stress
    { targets: { S: `require` }, model_targets: `EF_G`, expected: false },
    { targets: { S: `require` }, model_targets: `EFSH_G`, expected: true },
    // direct-only force/stress excludes gradient-based and force-less models
    { fs_mode: `direct`, model_targets: `EFS_G`, expected: false },
    { fs_mode: `direct`, model_targets: `EFS_DM`, expected: true },
    { fs_mode: `gradient`, targets: {}, model_targets: `E`, expected: false },
  ] as const)(
    `targets filter $targets/$fs_mode matches $model_targets -> $expected`,
    ({ targets, fs_mode, model_targets, expected }) => {
      const filters = make_table_filters()
      if (targets) filters.targets = { ...targets }
      if (fs_mode) filters.fs_mode = fs_mode
      expect(filters.matches({ training_set: [], targets: model_targets })).toBe(expected)
    },
  )

  it(`round-trips the targets URL param incl. the cleared (empty) state`, () => {
    const filters = make_table_filters()
    filters.read(new URLSearchParams(`targets=-M,direct`))
    expect(filters.targets).toStrictEqual({ M: `exclude` })
    expect(filters.fs_mode).toBe(`direct`)
    expect(filters.targets_param).toBe(`-M,direct`)

    // absent param = default (require F); present-but-empty = no targets filter
    filters.read(new URLSearchParams(``))
    expect(filters.targets).toStrictEqual({ F: `require` })
    filters.read(new URLSearchParams(`targets=`))
    expect(filters.targets).toStrictEqual({})
    expect(filters.targets_param).toBe(``)
  })

  it(`set_training cycles require -> off and toggle_openness keeps at least one`, () => {
    const filters = make_table_filters()
    filters.set_training(`MPtrj`, `require`)
    expect(filters.training.MPtrj).toBe(`require`)
    filters.set_training(`MPtrj`, `exclude`)
    expect(filters.training.MPtrj).toBe(`exclude`)
    filters.set_training(`MPtrj`, `exclude`) // same mode again clears it
    expect(filters.training.MPtrj).toBeUndefined()

    for (const openness of [`OSOD`, `OSCD`, `CSOD`, `CSCD`] as const) {
      filters.toggle_openness(openness)
    }
    expect(filters.openness).toHaveLength(1) // refuses to hide the last value
    filters.clear()
    expect(filters.n_active).toBe(0)
  })
})

describe(`find_best_model`, () => {
  const make_model = (
    model_name: string,
    f1: unknown,
    overrides: Record<string, unknown> = {},
  ) =>
    ({
      model_name,
      training_set: [`MPtrj`],
      openness: `OSOD`,
      metrics: { discovery: { full_test_set: { F1: f1 } } },
      ...overrides,
    }) as unknown as ModelData

  it(`picks the model with the highest full-test-set F1`, () => {
    const models = [
      make_model(`low`, 0.5),
      make_model(`high`, 0.9),
      make_model(`mid`, 0.7),
    ]
    expect(find_best_model(models)?.model_name).toBe(`high`)
  })

  it(`returns null instead of a truthy empty object when no model qualifies`, () => {
    // regression: best_model used to seed with a truthy {}, so when no model qualified
    // the empty object rendered 'undefined' model names instead of being falsy
    expect(find_best_model([])).toBeNull()
  })

  it(`ranks by the requested discovery_set`, () => {
    const make_discovery_model = (model_name: string, full: number, uniq: number) =>
      make_model(model_name, full, {
        metrics: {
          discovery: { full_test_set: { F1: full }, unique_prototypes: { F1: uniq } },
        },
      })
    const models = [
      make_discovery_model(`best-full`, 0.9, 0.4),
      make_discovery_model(`best-uniq`, 0.5, 0.8),
    ]
    expect(find_best_model(models)?.model_name).toBe(`best-full`)
    expect(find_best_model(models, `unique_prototypes`)?.model_name).toBe(`best-uniq`)
  })

  it(`skips models with missing or non-numeric F1`, () => {
    const models = [
      make_model(`no-f1`, undefined),
      make_model(`nan-f1`, Number.NaN),
      make_model(`string-f1`, `0.99`),
      make_model(`no-discovery`, 0, { metrics: {} }),
      make_model(`valid`, 0.6),
    ]
    expect(find_best_model(models)?.model_name).toBe(`valid`)
    expect(find_best_model(models.slice(0, 4))).toBeNull()
  })

  it(`finds a best model in the real MODELS data`, () => {
    const best = find_best_model(MODELS)
    if (!best) throw new Error(`expected a best model in real data`)
    const discovery = best.metrics?.discovery
    const f1 = typeof discovery === `object` ? discovery.full_test_set?.F1 : undefined
    expect(f1).toBeGreaterThan(0.5)
  })
})

describe(`ALL_TRAINING_SETS`, () => {
  it(`lists only datasets used by at least one model, in datasets.yml order`, () => {
    expect(ALL_TRAINING_SETS.length).toBeGreaterThan(5)
    for (const key of [`MPtrj`, `OMat24`, `sAlex`, `MP 2022`]) {
      expect(ALL_TRAINING_SETS, `missing ${key}`).toContain(key)
    }
    const used_keys = new Set<string>(MODELS.flatMap((model) => model.training_set))
    expect(ALL_TRAINING_SETS.every((key) => used_keys.has(key))).toBe(true)
  })
})

// NB: CPS_CONFIG defaults + reactivity are covered in combined-scores.test.ts
describe(`update_models_cps`, () => {
  it(`should update CPS for models based on metrics and current weights`, () => {
    // Skip test if no models available
    if (MODELS.length === 0) return

    // Act: Call the function under test
    update_models_cps(MODELS, CPS_CONFIG)

    // Assert: Check if at least one model has a non-NaN CPS value
    const models_with_cps = MODELS.filter((model) => !isNaN(Number(model.CPS)))
    expect(models_with_cps.length).toBeGreaterThan(0)
  })

  it(`should set CPS to NaN when required metrics are missing`, () => {
    // Skip test if no models available
    if (MODELS.length === 0) return

    // Act: Set weights that would require all metrics to be present
    CPS_CONFIG.F1.weight = 0.3
    CPS_CONFIG.RMSD.weight = 0.3
    CPS_CONFIG.κ_SRME.weight = 0.4
    update_models_cps(MODELS, CPS_CONFIG)

    // Assert: Check if some models have NaN CPS (due to missing metrics)
    const models_with_nan_cps = MODELS.filter((model) => isNaN(Number(model.CPS)))
    expect(models_with_nan_cps.length).toBeGreaterThan(0)
  })

  it(`computes CPS from the same RMSD symprec the table displays`, () => {
    // regression: update_models_cps read symprec=1e-5 RMSD while the table shows
    // symprec=1e-2 (ALL_METRICS.RMSD.path), so CPS could silently diverge from the
    // displayed value whenever the two symprec levels' RMSDs differ
    expect(ALL_METRICS.RMSD.path).toBe(`metrics.geo_opt.symprec=1e-2`)

    const model = {
      metrics: {
        discovery: { unique_prototypes: { F1: 0.8 } },
        geo_opt: { 'symprec=1e-2': { rmsd: 0.01 }, 'symprec=1e-5': { rmsd: 0.2 } },
        phonons: { kappa_103: { κ_SRME: 0.5 } },
      },
    } as unknown as ModelData
    const rmsd_only_config = {
      ...DEFAULT_CPS_CONFIG,
      F1: { ...DEFAULT_CPS_CONFIG.F1, weight: 0 },
      RMSD: { ...DEFAULT_CPS_CONFIG.RMSD, weight: 1 },
      κ_SRME: { ...DEFAULT_CPS_CONFIG.κ_SRME, weight: 0 },
    }
    update_models_cps([model], rmsd_only_config)
    expect(model.CPS).toBe(calculate_cps(undefined, 0.01, undefined, rmsd_only_config))
    expect(model.CPS).not.toBe(calculate_cps(undefined, 0.2, undefined, rmsd_only_config))
  })

  it(`should handle different weight configurations correctly`, () => {
    // First configuration - only F1 matters
    CPS_CONFIG.F1.weight = 1.0
    CPS_CONFIG.RMSD.weight = 0.0
    CPS_CONFIG.κ_SRME.weight = 0.0
    update_models_cps(MODELS, CPS_CONFIG)

    const f1_cps_values = MODELS.map((model) => Number(model.CPS))

    // Second configuration - only RMSD matters
    CPS_CONFIG.F1.weight = 0.0
    CPS_CONFIG.RMSD.weight = 1.0
    CPS_CONFIG.κ_SRME.weight = 0.0
    update_models_cps(MODELS, CPS_CONFIG)

    const rmsd_cps_values = MODELS.map((model) => Number(model.CPS))

    expect(f1_cps_values).not.toStrictEqual(rmsd_cps_values)

    CPS_CONFIG.F1.weight = 0.5
    CPS_CONFIG.RMSD.weight = 0
    CPS_CONFIG.κ_SRME.weight = 0.5
    update_models_cps(MODELS, CPS_CONFIG)

    const f1_rmsd_cps_values = MODELS.map((model) => Number(model.CPS))
    expect(f1_rmsd_cps_values).not.toStrictEqual(f1_cps_values)
    expect(f1_rmsd_cps_values).not.toStrictEqual(rmsd_cps_values)
  })
})

describe(`fig-helpers payload styling`, () => {
  it.each([
    [`ascending key`, (mdl: { v: number }) => mdl.v, [1, 2, 3]],
    [`descending via negation`, (mdl: { v: number }) => -mdl.v, [3, 2, 1]],
  ])(
    `order_models sorts by %s without mutating the input`,
    (_label, key_fn, expected) => {
      const input = [{ v: 2 }, { v: 3 }, { v: 1 }]
      expect(order_models(input, key_fn).map((mdl) => mdl.v)).toEqual(expected)
      expect(input.map((mdl) => mdl.v)).toEqual([2, 3, 1]) // input untouched
    },
  )

  it(`attach_style attaches MODELS colors and sorts models by discovery F1 desc`, () => {
    const keys = [`eSEN-30m-oam`, `equiformer-v3-oam`, `chgnet-0.3.0`]
    const f1 = (key: string) => {
      const disc = MODELS.find((mdl) => mdl.model_key === key)?.metrics?.discovery
      return typeof disc === `object`
        ? (disc?.unique_prototypes?.F1 ?? -Infinity)
        : -Infinity
    }
    const styled = attach_style({ shared: 1, models: keys.map((key) => ({ key })) })

    expect(styled.shared).toBe(1) // non-model shared fields are preserved
    const ordered_f1 = styled.models.map((mdl) => f1(mdl.key))
    expect(ordered_f1).toEqual([...ordered_f1].sort((row_a, row_b) => row_b - row_a))
    for (const mdl of styled.models) {
      expect(mdl.color).toBe(MODELS.find((model) => model.model_key === mdl.key)?.color)
    }
  })
})

// js-yaml >= 4.2 (YAML 1.2 core schema) parses underscore-grouped integers like
// `154_719` as strings instead of numbers, silently breaking numeric consumers
// (dataset sizes, model_params, ...). Parse each model/data YAML and flag any
// scalar that stayed a string only because it looks like an underscore-grouped
// integer, catching the regression that hit CI on PR #331 in any context (mapping
// value, sequence item or flow collection). Use plain integers (`154719`) instead.
const repo_root = path.resolve(import.meta.dirname, `..`, `..`, `..`)

const yaml_files = [`models`, `data`].flatMap((dir) =>
  readdirSync(path.join(repo_root, dir), { recursive: true })
    .map(String)
    .filter((name) => name.endsWith(`.yml`))
    .map((name) => path.join(dir, name)),
)

// a value js-yaml keeps as a string solely because YAML 1.2 forbids `_` in ints
const underscore_int = /^[-+]?\d{1,3}(?:_\d{3})+$/

function find_underscore_numbers(node: unknown, node_path: string): string[] {
  if (typeof node === `string`) {
    return underscore_int.test(node) ? [`${node_path} = ${node}`] : []
  }
  if (Array.isArray(node)) {
    return node.flatMap((item, idx) =>
      find_underscore_numbers(item, `${node_path}[${idx}]`),
    )
  }
  if (node && typeof node === `object`) {
    return Object.entries(node).flatMap(([key, value]) =>
      find_underscore_numbers(value, node_path ? `${node_path}.${key}` : key),
    )
  }
  return []
}

describe(`YAML data files use plain integers (no underscore separators)`, () => {
  it.each(yaml_files)(`%s`, (rel_path) => {
    const data = yaml_load(readFileSync(path.join(repo_root, rel_path), `utf-8`))
    const offenders = find_underscore_numbers(data, ``)
    expect(offenders, `${rel_path} has underscore-separated numbers`).toEqual([])
  })
})
