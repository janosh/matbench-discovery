import type { TargetType } from '$lib'
import { model_is_compliant } from '$lib'
import {
  create_combined_filter,
  format_date,
  format_train_set,
  get_geo_opt_property,
  targets_tooltips,
} from '$lib/metrics-table-helpers'
import type { ModelData } from '$lib/types'
import { beforeEach, describe, expect, it, vi } from 'vitest'

// Mock the imports that the helper functions depend on
vi.mock(`$lib`, () => ({
  model_is_compliant: vi.fn(),
  get_pred_file_urls: vi.fn().mockReturnValue([]),
}))

vi.mock(`$data/datasets.yml`, () => ({
  default: {
    'MP 2022': {
      title: `Materials Project 2022`,
      url: `https://materialsproject.org`,
      n_structures: 100000,
      n_materials: 50000,
    },
    MPtrj: {
      title: `Materials Project Trajectories`,
      url: `https://materialsproject.org/trajectories`,
      n_structures: 200000,
      n_materials: 75000,
    },
    'Custom Set': {
      title: `Custom Dataset`,
      url: `https://example.com`,
      n_structures: 10000,
    },
  },
}))

describe(`metrics-table-helpers`, () => {
  describe(`targets_tooltips`, () => {
    it.each([
      [`E`, `Energy`],
      [`EF_G`, `Energy with gradient-based forces`],
      [`EFS_DM`, `Energy with direct forces, stress, and magmoms`],
      [`EF_D`, `Energy with direct forces`],
      [`EFS_G`, `Energy with gradient-based forces and stress`],
    ])(`contains tooltip for %s target type`, (target, expected) => {
      expect(targets_tooltips[target as TargetType]).toBe(expected)
    })

    it(`contains all expected tooltip keys`, () => {
      expect(Object.keys(targets_tooltips).length).toBe(7)
    })
  })

  describe(`format_long_date`, () => {
    it(`formats date in long format`, () => {
      // Use a fixed date for testing
      const date = `2023-05-15`

      // Mock Date to return consistent results
      const original_date = Date
      const mock_date = class extends Date {
        constructor(date_str?: string | number | Date) {
          if (date_str) {
            super(date_str as string | number)
          } else {
            super(`2023-05-15T12:00:00Z`) // Mock current date
          }
        }

        toLocaleDateString(): string {
          return `Monday, May 15, 2023`
        }
      }

      // Override Date constructor
      Object.defineProperty(globalThis, `Date`, {
        value: mock_date,
        writable: true,
      })

      expect(format_date(date)).toBe(`Monday, May 15, 2023`)

      // Restore original Date
      Object.defineProperty(globalThis, `Date`, {
        value: original_date,
        writable: true,
      })
    })
  })

  describe(`format_train_set`, () => {
    it.each([
      {
        case: `single training set`,
        input: [`MP 2022`],
        expected_contains: [
          `<span title="`,
          `data-sort-value="50000"`,
          `<a href="https://materialsproject.org" target="_blank" rel="noopener noreferrer">MP 2022</a>`,
          `50k`,
          `50,000 materials in training set`,
        ],
      },
      {
        case: `multiple training sets`,
        input: [`MP 2022`, `MPtrj`],
        expected_contains: [
          `<a href="https://materialsproject.org" target="_blank" rel="noopener noreferrer">MP 2022</a>`,
          `<a href="https://materialsproject.org/trajectories" target="_blank" rel="noopener noreferrer">MPtrj</a>`,
          `MP 2022</a>+<a`,
          `data-sort-value="125000"`,
          `125k`,
          `&#013;• Materials Project 2022: 50,000 materials`,
          `• Materials Project Trajectories: 75,000 materials`,
        ],
      },
      {
        case: `training set with materials and structures`,
        input: [`MP 2022`],
        expected_contains: [
          `50k <small>(100k)</small>`,
          `50,000 materials in training set (100,000 structures`,
        ],
      },
      {
        case: `training set without n_materials`,
        input: [`Custom Set`],
        expected_contains: [
          `data-sort-value="10000"`,
          `10k`,
          `<a href="https://example.com" target="_blank" rel="noopener noreferrer">Custom Set</a>`,
        ],
        not_contains: [`<small>`],
      },
    ])(`formats $case correctly`, ({ input, expected_contains, not_contains }) => {
      const result = format_train_set(input)

      // Check for expected content
      for (const content of expected_contains) {
        expect(result).toContain(content)
      }

      // Check for content that should not be present
      if (not_contains) {
        for (const content of not_contains) {
          expect(result).not.toContain(content)
        }
      }
    })

    it(`handles missing training sets gracefully with warnings`, () => {
      // Mock console.warn
      const console_spy = vi.spyOn(console, `warn`).mockImplementation(() => {})

      const result = format_train_set([`MP 2022`, `NonExistent`])

      // Should warn about missing training set with exact message
      expect(console_spy).toHaveBeenCalledWith(
        `Training set NonExistent not found in TRAINING_SETS`,
      )

      // Should still format the existing training set correctly
      expect(result).toContain(
        `<a href="https://materialsproject.org" target="_blank" rel="noopener noreferrer">MP 2022</a>`,
      )
      expect(result).toContain(`50k`)

      // Should not include the missing dataset name anywhere
      expect(result).not.toContain(`NonExistent`)

      console_spy.mockRestore()
    })

    it(`shows both materials and structures when they differ`, () => {
      const result = format_train_set([`MP 2022`])

      // Verify both material and structure counts are shown
      expect(result).toContain(`50k <small>(100k)</small>`)

      // Check tooltip includes both counts with proper formatting
      expect(result).toContain(
        `50,000 materials in training set (100,000 structures counting all DFT relaxation`,
      )
    })

    it(`formats training sets without n_materials correctly using n_structures`, () => {
      const result = format_train_set([`Custom Set`])

      // Should use n_structures as n_materials
      expect(result).toContain(`data-sort-value="10000"`)
      expect(result).toContain(`10k`)

      // Should not contain small tag since materials and structures are the same
      expect(result).not.toContain(`<small>`)

      // Should include Custom Dataset in the result
      expect(result).toContain(
        `<a href="https://example.com" target="_blank" rel="noopener noreferrer">Custom Set</a>`,
      )
    })
  })

  describe(`get_geo_opt_property`, () => {
    const geo_opt = {
      'symprec=0.1': {
        rmsd: 0.025,
        energy_diff: 0.01,
      },
    }

    it.each([
      {
        case: `valid property`,
        geo_opt,
        symprec: `0.1`,
        property: `rmsd`,
        expected: 0.025,
      },
      {
        case: `another valid property`,
        geo_opt,
        symprec: `0.1`,
        property: `energy_diff`,
        expected: 0.01,
      },
      {
        case: `missing symprec`,
        geo_opt,
        symprec: `0.2`,
        property: `rmsd`,
        expected: undefined,
      },
      {
        case: `missing property`,
        geo_opt,
        symprec: `0.1`,
        property: `missing_prop`,
        expected: undefined,
      },
    ])(`returns $expected for $case`, ({ geo_opt, symprec, property, expected }) => {
      const result = get_geo_opt_property<number>(geo_opt, symprec, property)
      expect(result).toBe(expected)
    })

    it.each([
      [`null input`, null],
      [`string input`, `not an object`],
      [`undefined input`, undefined],
    ])(`handles %s gracefully`, (_case, input) => {
      const result = get_geo_opt_property<number>(input, `0.1`, `rmsd`)
      expect(result).toBeUndefined()
    })
  })

  describe(`create_combined_filter`, () => {
    // Reset the mock for each test
    beforeEach(() => {
      vi.resetAllMocks()
    })

    it.each([
      {
        case: `user filter returns false`,
        model_filter_returns: false,
        show_energy: true,
        show_noncomp: true,
        model: { targets: `E` },
        expected_result: false,
        should_check_compliance: false,
      },
      {
        case: `energy model with show_energy=false`,
        model_filter_returns: true,
        show_energy: false,
        show_noncomp: true,
        model: { targets: `E` },
        expected_result: false,
        should_check_compliance: false,
      },
      {
        case: `energy model with show_energy=true`,
        model_filter_returns: true,
        show_energy: true,
        show_noncomp: true,
        model: { targets: `E` },
        expected_result: true,
        should_check_compliance: true,
        is_compliant: true,
      },
      {
        case: `force model`,
        model_filter_returns: true,
        show_energy: false,
        show_noncomp: true,
        model: { targets: `EF_G` },
        expected_result: true,
        should_check_compliance: true,
        is_compliant: true,
      },
      {
        case: `non-compliant model with show_noncomp=false`,
        model_filter_returns: true,
        show_energy: true,
        show_noncomp: false,
        model: { targets: `EF_G` },
        expected_result: false,
        should_check_compliance: true,
        is_compliant: false,
      },
      {
        case: `non-compliant model with show_noncomp=true`,
        model_filter_returns: true,
        show_energy: true,
        show_noncomp: true,
        model: { targets: `EF_G` },
        expected_result: true,
        should_check_compliance: true,
        is_compliant: false,
      },
      {
        case: `all conditions pass`,
        model_filter_returns: true,
        show_energy: true,
        show_noncomp: true,
        model: { targets: `EF_G` },
        expected_result: true,
        should_check_compliance: true,
        is_compliant: true,
      },
    ])(
      `$case`,
      ({
        model_filter_returns,
        show_energy,
        show_noncomp,
        model,
        expected_result,
        should_check_compliance,
        is_compliant,
      }) => {
        const mock_model_filter = vi.fn().mockReturnValue(model_filter_returns)

        if (should_check_compliance && is_compliant !== undefined) {
          vi.mocked(model_is_compliant).mockReturnValue(is_compliant)
        }

        const filter = create_combined_filter(
          mock_model_filter,
          show_energy,
          show_noncomp,
        )

        expect(filter(model as ModelData)).toBe(expected_result)
        expect(mock_model_filter).toHaveBeenCalledWith(model)

        if (should_check_compliance) {
          expect(model_is_compliant).toHaveBeenCalledWith(model)
        } else {
          expect(model_is_compliant).not.toHaveBeenCalled()
        }
      },
    )
  })
})
