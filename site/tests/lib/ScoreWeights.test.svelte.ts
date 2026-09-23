import ScoreWeights from '$lib/ScoreWeights.svelte'
import {
  DEFAULT_CDS_CONFIG,
  DEFAULT_CMDS_CONFIG,
  DEFAULT_CPS_CONFIG,
  type CpsConfig,
} from '$lib/combined-scores.svelte'
import { ALL_METRICS, MD_METRICS } from '$lib/labels'
import { flushSync } from 'svelte'
import { describe, expect, it } from 'vite-plus/test'
import { doc_query, mount } from '../index'

describe(`ScoreWeights`, () => {
  it.each([
    [DEFAULT_CPS_CONFIG, ALL_METRICS.CPS, [`50`, `40`, `10`]],
    [DEFAULT_CMDS_CONFIG, MD_METRICS.md_combined_score, [`30`, `20`, `20`, `30`]],
    [
      DEFAULT_CDS_CONFIG,
      ALL_METRICS.diatomics_combined_score,
      [`44.4`, `22.2`, `11.1`, `22.2`],
    ],
  ] as const)(
    `edits, normalizes and resets configuration %#`,
    (default_config, title_label, percentages) => {
      const config = $state(structuredClone(default_config))
      mount(ScoreWeights, {
        target: document.body,
        props: { config, default_config, title_label },
      })
      const inputs = [
        ...document.querySelectorAll<HTMLInputElement>(`.score-weights input`),
      ]
      const reset = doc_query<HTMLButtonElement>(`.score-weights button`)
      expect(inputs.map(({ value }) => value)).toEqual(percentages)
      expect(
        inputs.every((input) => input.getAttribute(`aria-label`)?.endsWith(`weight (%)`)),
      ).toBe(true)
      expect(reset.disabled).toBe(true)

      const edit = (value: string) => {
        inputs[0].value = value
        inputs[0].dispatchEvent(new Event(`change`, { bubbles: true }))
        flushSync()
        return Object.values(config).map(({ weight }) => weight)
      }
      const original = Object.values(default_config).map(({ weight }) => weight)
      const edited = edit(`75`)
      expect(edited[0]).toBe(0.75)
      for (let idx = 1; idx < edited.length; idx++)
        expect(edited[idx]).toBeCloseTo((0.25 * original[idx]) / (1 - original[0]), 12)
      expect(edited.reduce((sum, weight) => sum + weight, 0)).toBeCloseTo(1, 12)
      expect(reset.disabled).toBe(false)
      expect(edit(`100`)).toEqual([1, ...original.slice(1).map(() => 0)])
      // Leaving a single-component score redistributes the remainder equally.
      const spread = edit(`25`)
      expect(spread[0]).toBe(0.25)
      for (const weight of spread.slice(1))
        expect(weight).toBeCloseTo(0.75 / (spread.length - 1), 12)
      expect(edit(`0`)[0]).toBe(0)

      reset.click()
      flushSync()
      expect(Object.values(config).map(({ weight }) => weight)).toEqual(original)
      expect(inputs.map(({ value }) => value)).toEqual(percentages)
      expect(reset.disabled).toBe(true)
    },
  )

  it.each([``, `-1`, `101`, `1e309`])(`restores invalid percentage %s`, (value) => {
    const config = $state<CpsConfig>(structuredClone(DEFAULT_CPS_CONFIG))
    mount(ScoreWeights, { target: document.body, props: { config } })
    const input = doc_query<HTMLInputElement>(`.score-weights input`)
    input.value = value
    input.dispatchEvent(new Event(`change`, { bubbles: true }))
    flushSync()
    expect(input.value).toBe(`50`)
    expect(config).toEqual(DEFAULT_CPS_CONFIG)
    // External URL restoration updates the same controls without local draft state.
    config.F1.weight = 0.7
    config.κ_SRME.weight = 0.2
    flushSync()
    expect(input.value).toBe(`70`)
  })
})
