import { default as DATASETS } from '$data/datasets.yml'
import { ALL_TRAINING_SETS } from './models.svelte'
import type { FilterPreset } from './url-state.svelte'

// Built-in presets shipped with the site. `Compliant` recreates the pre-2026 leaderboard
// "compliant models" cohort in one click: open source + open data (OSOD) models trained
// exclusively on MP-anchored datasets (those with compliant: true in datasets.yml, same
// source of truth as the Python API's model_is_compliant).
export const BUILTIN_PRESETS: Record<string, FilterPreset> = {
  Compliant: {
    training: Object.fromEntries(
      ALL_TRAINING_SETS.filter((key) => !DATASETS[key]?.compliant).map((key) => [
        key,
        `exclude` as const,
      ]),
    ),
    openness: [`OSOD`],
    description:
      `Open source + open data models trained exclusively on MP-anchored datasets ` +
      `(the former "compliant" leaderboard cohort)`,
  },
}

const STORAGE_KEY = `metrics-table-filter-presets`

function load_user_presets(): Record<string, FilterPreset> {
  try {
    // globalThis guard: localStorage doesn't exist during SSR/prerendering
    return JSON.parse(globalThis.localStorage?.getItem(STORAGE_KEY) ?? `{}`)
  } catch {
    return {}
  }
}

// User-defined presets, persisted to localStorage and shared by all metrics tables
export const user_presets = $state<Record<string, FilterPreset>>(load_user_presets())

function persist(): void {
  // setItem can throw (Safari private mode, quota exceeded); losing persistence
  // shouldn't crash the save/delete click handlers
  try {
    globalThis.localStorage?.setItem(STORAGE_KEY, JSON.stringify(user_presets))
  } catch (error) {
    console.error(`Failed to persist filter presets:`, error)
  }
}

export function save_user_preset(name: string, preset: FilterPreset): void {
  user_presets[name] = preset
  persist()
}

export function delete_user_preset(name: string): void {
  Reflect.deleteProperty(user_presets, name)
  persist()
}
