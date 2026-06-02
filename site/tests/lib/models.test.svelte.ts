import { CPS_CONFIG } from '$lib/combined_perf_score.svelte'
import {
  calculate_training_sizes,
  COMPLIANT_TRAINING_SETS,
  model_is_compliant,
  MODEL_METADATA_PATHS,
  MODELS,
  update_models_cps,
} from '$lib/models.svelte'
import per_elem_each_errors from '$routes/models/per-element-each-errors.json'
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
  it(`should be defined and be an array`, () => {
    expect(MODELS).toBeDefined()
    expect(Array.isArray(MODELS)).toBe(true)
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
    expect(MODEL_METADATA_PATHS).toBeDefined()
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

describe(`model_is_compliant`, () => {
  it.each([
    { training_set: [`MPtrj`], openness: `OSOD`, expected: true },
    { training_set: [`MP 2022`], openness: `OSOD`, expected: true },
    { training_set: [`MPF`], openness: `OSOD`, expected: true },
    { training_set: [`MP Graphs`], openness: `OSOD`, expected: true },
    { training_set: [`MPtrj`, `MP 2022`], openness: `OSOD`, expected: true },
    // MDR PBE Phonons in MPtrj is compliant (phonon fine-tuning on MPtrj subset)
    { training_set: [`MDR PBE Phonons in MPtrj`], openness: `OSOD`, expected: true },
    {
      training_set: [`MPtrj`, `MDR PBE Phonons in MPtrj`],
      openness: `OSOD`,
      expected: true,
    },
  ])(
    `returns true for compliant training_set=$training_set`,
    ({ training_set, openness, expected }) => {
      const model = { training_set, openness } as Parameters<typeof model_is_compliant>[0]
      expect(model_is_compliant(model)).toBe(expected)
    },
  )

  it.each([
    { training_set: [`OMat24`], openness: `OSOD`, expected: false },
    { training_set: [`GNoME`], openness: `OSOD`, expected: false },
    { training_set: [`MPtrj`, `OMat24`], openness: `OSOD`, expected: false },
    { training_set: [`Alex`], openness: `OSOD`, expected: false },
  ])(
    `returns false for non-compliant training_set=$training_set`,
    ({ training_set, openness, expected }) => {
      const model = { training_set, openness } as Parameters<typeof model_is_compliant>[0]
      expect(model_is_compliant(model)).toBe(expected)
    },
  )

  it.each([
    { training_set: [`MPtrj`], openness: `CSOD`, expected: false },
    { training_set: [`MPtrj`], openness: `OSCD`, expected: false },
    { training_set: [`MPtrj`], openness: `CSCD`, expected: false },
  ])(
    `returns false for non-OSOD openness=$openness`,
    ({ training_set, openness, expected }) => {
      const model = { training_set, openness } as Parameters<typeof model_is_compliant>[0]
      expect(model_is_compliant(model)).toBe(expected)
    },
  )

  it(`defaults to OSOD when openness is undefined`, () => {
    const model = { training_set: [`MPtrj`] } as Parameters<typeof model_is_compliant>[0]
    expect(model_is_compliant(model)).toBe(true)
  })
})

describe(`COMPLIANT_TRAINING_SETS`, () => {
  it(`matches expected compliant datasets from datasets.yml`, () => {
    // This test ensures Python and TypeScript compute the same compliant sets
    // Both should derive from datasets.yml where compliant: true
    const expected = [`MP 2022`, `MPtrj`, `MPF`, `MP Graphs`, `MDR PBE Phonons in MPtrj`]
    expect(COMPLIANT_TRAINING_SETS.toSorted()).toStrictEqual(expected.toSorted())
  })

  it(`is exported as an array`, () => {
    expect(Array.isArray(COMPLIANT_TRAINING_SETS)).toBe(true)
  })
})

describe(`CPS_CONFIG`, () => {
  it(`should be defined and match DEFAULT_CPS_CONFIG initially`, () => {
    expect(CPS_CONFIG).toBeDefined()
    expect(CPS_CONFIG.F1.weight).toBeDefined()
    expect(CPS_CONFIG.RMSD.weight).toBeDefined()
    expect(CPS_CONFIG.κ_SRME.weight).toBeDefined()
  })

  it(`should be reactive (modifiable)`, () => {
    // Store original weights
    const original_f1_weight = CPS_CONFIG.F1.weight

    // Modify weight
    CPS_CONFIG.F1.weight = 0.8

    // Check that it was updated
    expect(CPS_CONFIG.F1.weight).toBe(0.8)

    // Restore original weight
    CPS_CONFIG.F1.weight = original_f1_weight
  })
})

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
