import { afterNavigate, replaceState } from '$app/navigation'
import { page } from '$app/state'
import * as d3_sc from 'd3-scale-chromatic'
import type { D3InterpolateName } from 'matterviz/colors'
import type { SortDir } from './types'

type PageState = Parameters<typeof replaceState>[1]
export type SortState = { column: string; dir: SortDir }
export type UrlParamEntry = [key: string, value: string, default_value?: string]
type ValidValues<T extends string> = ReadonlySet<T> | Record<string, unknown>
const sort_dirs = new Set<SortDir>([`asc`, `desc`])

// Boolean flag params: default-true flags encode "off" as `0` (any other value = on),
// default-false flags encode "on" as `1` (any other value = off). At their default,
// flags are omitted from the URL.
export const bool_from_param = (
  params: URLSearchParams,
  key: string,
  fallback = false,
): boolean => (fallback ? params.get(key) !== `0` : params.get(key) === `1`)

export const bool_url_entry = (
  key: string,
  value: boolean,
  fallback = false,
): UrlParamEntry => [key, value === fallback ? `` : value ? `1` : `0`]

export function valid_query_param<T extends string>(
  params: URLSearchParams,
  key: string,
  fallback: T,
  valid_values: ValidValues<T>,
): T {
  const value = params.get(key)
  if (!value) return fallback

  const is_valid =
    valid_values instanceof Set
      ? valid_values.has(value)
      : Object.hasOwn(valid_values, value)
  return is_valid ? (value as T) : fallback
}

// valid_columns is optional: the sortable-column set lives inside the table component
// (unknown columns are a harmless no-op there); pass it where known to reject garbage
export const sort_from_query = (
  params: URLSearchParams,
  default_sort: SortState,
  valid_columns?: ValidValues<string>,
): SortState => ({
  column: valid_columns
    ? valid_query_param(params, `sort`, default_sort.column, valid_columns)
    : (params.get(`sort`) ?? default_sort.column),
  dir: valid_query_param(params, `dir`, default_sort.dir, sort_dirs),
})

// write-side counterpart of sort_from_query
export const sort_url_entries = (
  sort: SortState,
  default_sort: SortState,
): UrlParamEntry[] => [
  [`sort`, sort.column, default_sort.column],
  [`dir`, sort.dir, default_sort.dir],
]

// -- Weighted-score radar weights as a single URL param ------------------------
// Serialized as comma-joined values in config-key order, e.g. weights=0.5,0.4,0.1.
type WeightsConfig = Record<string, { weight: number }>

const round_weight = (weight: number): number => Math.round(weight * 1000) / 1000

// Empty string when weights match the defaults (so sync_url_params drops the param)
export function weights_to_param(
  config: WeightsConfig,
  default_config: WeightsConfig,
): string {
  const keys = Object.keys(config)
  const is_default = keys.every(
    (key) =>
      round_weight(config[key].weight) ===
      round_weight(default_config[key]?.weight ?? NaN),
  )
  return is_default ? `` : keys.map((key) => round_weight(config[key].weight)).join(`,`)
}

// Parse a weights param and write it into config (normalized to sum 1). A missing
// OR malformed param (wrong count, negative/non-finite, all-zero) resets to
// default_config: weight configs are shared module state, so without the reset,
// weights customized earlier in the session would survive navigating to a
// weights-less (or mangled) URL while all other URL-bound page state (sort, axes)
// resets - and the URL-sync effect would then launder those stale weights back into
// a valid-looking URL.
export function apply_weights_param(
  param: string | null,
  config: WeightsConfig,
  default_config: WeightsConfig,
): void {
  const keys = Object.keys(config)
  if (param) {
    // empty segments parse to NaN (not Number(``) which is 0) so a mangled URL like
    // weights=0.5,,0.5 is rejected by the finiteness check instead of zeroing a metric
    const values = param.split(`,`).map((part) => (part.trim() ? Number(part) : NaN))
    const total = values.reduce((sum, val) => sum + val, 0)
    if (
      values.length === keys.length &&
      values.every((val) => Number.isFinite(val) && val >= 0) &&
      total > 0
    ) {
      for (const [idx, key] of keys.entries()) config[key].weight = values[idx] / total
      return
    }
  }
  for (const key of keys)
    config[key].weight = default_config[key]?.weight ?? config[key].weight
}

// color_scale param: valid d3 interpolate names, defaulting to Viridis
const d3_color_scale_names = new Set(
  Object.keys(d3_sc).filter((key) => key.startsWith(`interpolate`)),
) as Set<D3InterpolateName>

export const url_color_scale = {
  default: `interpolateViridis` as D3InterpolateName,
  read: (params: URLSearchParams): D3InterpolateName =>
    valid_query_param(
      params,
      `color_scale`,
      url_color_scale.default,
      d3_color_scale_names,
    ),
  entry: (value: D3InterpolateName): UrlParamEntry => [
    `color_scale`,
    value,
    url_color_scale.default,
  ],
}

// -- Metrics-table model filters (training data, openness, targets, heatmap) -------
// Encoded as:
//   train=MPtrj,-OMat24   comma list of dataset keys; a bare key requires the dataset
//                         in a model's training set (multiple keys AND together), a
//                         -prefixed key excludes models trained on that dataset.
//                         Omitted when no dataset is filtered.
//   openness=OSOD,OSCD    subset of openness values to show; omitted when all shown
//   targets=F,-M,direct   same require/exclude scheme for predicted outputs (F/S/M/H)
//                         plus an optional direct|gradient token restricting how F/S
//                         are computed. Omitted at the default (F required, i.e.
//                         energy-only models hidden); `targets=` (empty) = no filter.
//   heatmap=0             heatmap colors off (default on)
export const OPENNESS_OPTIONS = [`OSOD`, `OSCD`, `CSOD`, `CSCD`] as const
export type Openness = (typeof OPENNESS_OPTIONS)[number]
export const TRAIN_FILTER_MODES = [`require`, `exclude`] as const
export type TrainFilterMode = (typeof TRAIN_FILTER_MODES)[number]
// filterable predicted outputs (energy is universal, so not filterable): forces,
// stress, magmoms, Hessian. Keys match the letters in model.targets (e.g. EFS_GM).
export const TARGET_OUTPUTS = {
  F: `forces`,
  S: `stress`,
  M: `magmoms`,
  H: `Hessian`,
} as const
export type TargetOutput = keyof typeof TARGET_OUTPUTS
const target_output_keys = Object.keys(TARGET_OUTPUTS) as TargetOutput[]
// how forces/stress are computed: direct model heads vs energy gradients
export const FS_MODES = [`any`, `direct`, `gradient`] as const
export type FsMode = (typeof FS_MODES)[number]
const DEFAULT_TARGETS = { F: `require` } as const
export const DEFAULT_TARGETS_PARAM = `F`
// a saved filter combination (see $lib/filter-presets.svelte.ts)
export type FilterPreset = {
  training: Record<string, TrainFilterMode>
  openness: readonly Openness[]
  targets?: Partial<Record<TargetOutput, TrainFilterMode>> // absent = default (require F)
  fs_mode?: FsMode
  description?: string // tooltip, only set on built-in presets
}
// minimal structural model shape keeps this module decoupled from $lib/types
type FilterableModel = { training_set: string[]; openness?: string; targets?: string }
const is_one_of = <Value extends string>(
  options: readonly Value[],
  value: unknown,
): value is Value => typeof value === `string` && options.includes(value as Value)

// Split a model.targets string like `EFS_GM` into its predicted outputs and the
// force/stress computation mode: prefix letters are E/F/S/H outputs, the suffix
// holds G(radient)/D(irect) plus M when the model also predicts magmoms.
// Exported so filter UIs can tally models per output with the same semantics.
export function parse_targets(targets = ``): {
  outputs: Set<string>
  fs_mode: FsMode | null
} {
  const [prefix, suffix = ``] = targets.split(`_`)
  const outputs = new Set<string>(prefix)
  if (suffix.includes(`M`)) outputs.add(`M`)
  const fs_mode = suffix.includes(`D`)
    ? `direct`
    : suffix.includes(`G`)
      ? `gradient`
      : null
  return { outputs, fs_mode }
}

export class UrlTableFilters {
  // dataset key -> require/exclude; keys absent from the record are unfiltered
  training = $state<Record<string, TrainFilterMode>>({})
  openness = $state<Openness[]>([...OPENNESS_OPTIONS])
  // predicted-output constraints; forces required by default (hides energy-only models)
  targets = $state<Partial<Record<TargetOutput, TrainFilterMode>>>({
    ...DEFAULT_TARGETS,
  })
  fs_mode = $state<FsMode>(`any`)
  show_heatmap = $state(true)

  constructor(readonly training_sets: string[]) {}

  // number of active non-default constraints (drives filter-button badges)
  get n_active(): number {
    return (
      Object.keys(this.training).length +
      (this.openness.length < OPENNESS_OPTIONS.length ? 1 : 0) +
      (this.targets_param === DEFAULT_TARGETS_PARAM ? 0 : 1)
    )
  }

  matches = (model: FilterableModel): boolean => {
    if (!this.openness.includes((model.openness ?? `OSOD`) as Openness)) return false
    const { outputs, fs_mode } = parse_targets(model.targets)
    const outputs_ok = Object.entries(this.targets).every(
      ([key, mode]) => outputs.has(key) === (mode === `require`),
    )
    if (!outputs_ok) return false
    // direct/gradient also drops models without any force/stress prediction
    if (this.fs_mode !== `any` && fs_mode !== this.fs_mode) return false
    return Object.entries(this.training).every(
      ([key, mode]) => model.training_set.includes(key) === (mode === `require`),
    )
  }

  // toggle a dataset constraint; picking the already-active mode clears it
  set_training = (key: string, mode: TrainFilterMode): void => {
    if (this.training[key] === mode) Reflect.deleteProperty(this.training, key)
    else this.training[key] = mode
  }

  // toggle a predicted-output constraint, same cycling as set_training
  set_target = (key: TargetOutput, mode: TrainFilterMode): void => {
    if (this.targets[key] === mode) Reflect.deleteProperty(this.targets, key)
    else this.targets[key] = mode
  }

  // flip an openness value's membership (keeping canonical order), refusing to
  // hide the last one (would empty the table)
  toggle_openness = (value: Openness): void => {
    const next = OPENNESS_OPTIONS.filter((op) =>
      op === value ? !this.openness.includes(op) : this.openness.includes(op),
    )
    if (next.length > 0) this.openness = next
  }

  clear = (): void => {
    this.training = {}
    this.openness = [...OPENNESS_OPTIONS]
    this.targets = { ...DEFAULT_TARGETS }
    this.fs_mode = `any`
  }

  apply = (preset: FilterPreset): void => {
    // keep only known datasets + valid modes: stale localStorage presets (e.g. after a
    // dataset rename) would otherwise filter models invisibly - unrepresentable in the
    // URL (url_entries serializes canonical keys only) and not shown by any checkbox
    this.training = Object.fromEntries(
      Object.entries(preset.training).filter(
        ([key, mode]) =>
          this.training_sets.includes(key) && is_one_of(TRAIN_FILTER_MODES, mode),
      ),
    )
    // filter OPENNESS_OPTIONS (not spread the preset) to keep canonical order and
    // drop invalid tokens from hand-edited localStorage
    const shown = OPENNESS_OPTIONS.filter((op) => preset.openness.includes(op))
    this.openness = shown.length > 0 ? shown : [...OPENNESS_OPTIONS]
    this.targets = Object.fromEntries(
      Object.entries(preset.targets ?? DEFAULT_TARGETS).filter(
        ([key, mode]) =>
          is_one_of(target_output_keys, key) && is_one_of(TRAIN_FILTER_MODES, mode),
      ),
    )
    this.fs_mode = is_one_of(FS_MODES, preset.fs_mode) ? preset.fs_mode : `any`
  }

  // snapshot of the active filters, e.g. for saving as a preset
  get as_preset(): FilterPreset {
    return {
      training: { ...this.training },
      openness: [...this.openness],
      targets: { ...this.targets },
      fs_mode: this.fs_mode,
    }
  }

  read = (params: URLSearchParams): void => {
    const valid_sets = new Set(this.training_sets)
    const training: Record<string, TrainFilterMode> = {}
    for (const token of params.get(`train`)?.split(`,`) ?? []) {
      const exclude = token.startsWith(`-`)
      const key = exclude ? token.slice(1) : token
      if (valid_sets.has(key)) training[key] = exclude ? `exclude` : `require`
    }
    this.training = training

    const shown = params
      .get(`openness`)
      ?.split(`,`)
      .filter((token) => is_one_of(OPENNESS_OPTIONS, token))
    this.openness = shown?.length
      ? OPENNESS_OPTIONS.filter((op) => shown.includes(op))
      : [...OPENNESS_OPTIONS]

    // absent param = default (require F); present-but-empty `targets=` = no filter
    const targets_param = params.get(`targets`)
    const targets: Partial<Record<TargetOutput, TrainFilterMode>> = {}
    let fs_mode: FsMode = `any`
    for (const token of targets_param?.split(`,`).filter(Boolean) ?? []) {
      if (is_one_of(FS_MODES, token)) fs_mode = token
      else {
        const exclude = token.startsWith(`-`)
        const key = (exclude ? token.slice(1) : token) as TargetOutput
        if (key in TARGET_OUTPUTS) targets[key] = exclude ? `exclude` : `require`
      }
    }
    this.targets = targets_param === null ? { ...DEFAULT_TARGETS } : targets
    this.fs_mode = fs_mode

    this.show_heatmap = bool_from_param(params, `heatmap`, true)
  }

  // canonical serialization of the targets + fs_mode constraints (F,-M,direct)
  get targets_param(): string {
    return [
      ...target_output_keys
        .filter((key) => key in this.targets)
        .map((key) => (this.targets[key] === `exclude` ? `-${key}` : key)),
      ...(this.fs_mode === `any` ? [] : [this.fs_mode]),
    ].join(`,`)
  }

  get url_entries(): UrlParamEntry[] {
    // serialize in canonical dataset order so URLs are order-insensitive
    const train = this.training_sets
      .filter((key) => key in this.training)
      .map((key) => (this.training[key] === `exclude` ? `-${key}` : key))
      .join(`,`)
    const openness =
      this.openness.length < OPENNESS_OPTIONS.length ? this.openness.join(`,`) : ``
    return [
      [`train`, train],
      [`openness`, openness],
      [`targets`, this.targets_param, DEFAULT_TARGETS_PARAM],
      bool_url_entry(`heatmap`, this.show_heatmap, true),
    ]
  }
}

export function sync_url_params(entries: UrlParamEntry[], state: PageState): void {
  const params = new URLSearchParams(location.search)
  for (const [key, value, default_value = ``] of entries) {
    if (value === default_value) params.delete(key)
    else params.set(key, value)
  }

  // keep commas human-readable: they're legal in query strings (RFC 3986 sub-delims)
  // but URLSearchParams percent-encodes them, turning ?weights=0.5,0.4,0.1 into
  // ?weights=0.5%2C0.4%2C0.1. Decoding is value-preserving since URLSearchParams
  // parses literal and encoded commas identically
  const query = params.toString().replaceAll(`%2C`, `,`)
  const next_url = query ? `${location.pathname}?${query}` : location.pathname
  if (next_url !== `${location.pathname}${location.search}`) replaceState(next_url, state)
}

// Two-way URL query-param binding shared by all task pages. Reads state from the URL in
// afterNavigate (fires after the router is initialized, both on hydration and later
// navigations), then keeps the URL in sync with page state via replaceState. Gating
// writes on the first afterNavigate ensures the sync $effect never runs during the
// initial mount flush, which would throw "before router is initialized". Must be
// called during component init.
export function bind_url_params(
  read_params: ((params: URLSearchParams) => void) | null,
  entries: () => UrlParamEntry[],
): void {
  let url_ready = $state(false)

  afterNavigate(() => {
    read_params?.(page.url.searchParams)
    url_ready = true
  })

  $effect(() => {
    if (!url_ready) return
    sync_url_params(entries(), page.state)
  })
}
