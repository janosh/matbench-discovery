import DATASETS from '$data/datasets.yml'
import { storage_get_json, storage_set_json } from 'svelte-widgets/storage'
import { is_object } from 'svelte-widgets/utils'
import { ALL_TRAINING_SETS } from './models.svelte'
import type { FilterPreset } from './url-state.svelte'

// Built-in presets shipped with the site. `Compliant` recreates the pre-2026 leaderboard
// "compliant models" cohort in one click: open source + open data (OSOD) models trained
// exclusively on MP-anchored datasets (those with compliant: true in datasets.yml).
export const BUILTIN_PRESETS: Record<string, FilterPreset> = {
  Compliant: {
    training: Object.fromEntries(
      ALL_TRAINING_SETS.filter((key) => !DATASETS[key]?.compliant).map((key) => [
        key,
        `exclude` as const,
      ]),
    ),
    openness: [`OSOD`],
    description: `Open source + open data models trained exclusively on MP-anchored datasets (the former "compliant" leaderboard cohort)`,
  },
}

const STORAGE_KEY = `metrics-table-filter-presets`

const is_plain_object = (value: unknown): value is Record<string, unknown> =>
  is_object(value) && !Array.isArray(value)

// UrlTableFilters.apply validates modes and keys; reject only malformed shapes here.
const is_filter_preset = (value: unknown): value is FilterPreset =>
  is_plain_object(value) &&
  is_plain_object(value.training) &&
  Array.isArray(value.openness)

function load_user_presets(): Record<string, FilterPreset> {
  const stored = storage_get_json(STORAGE_KEY, {})
  if (!is_plain_object(stored)) return {}
  return Object.fromEntries(
    Object.entries(stored).filter((entry): entry is [string, FilterPreset] =>
      is_filter_preset(entry[1]),
    ),
  )
}

// User-defined presets, persisted to localStorage and shared by all metrics tables
export const user_presets = $state<Record<string, FilterPreset>>(load_user_presets())

const persist = (): void => storage_set_json(STORAGE_KEY, user_presets)

export function save_user_preset(name: string, preset: FilterPreset): void {
  user_presets[name] = preset
  persist()
}

export function delete_user_preset(name: string): void {
  Reflect.deleteProperty(user_presets, name)
  persist()
}
