import { afterNavigate, onNavigate, replaceState } from '$app/navigation'
import { page } from '$app/state'
import type { AfterNavigate } from '@sveltejs/kit'
import { untrack } from 'svelte'
import { is_d3_interpolate_name, type D3InterpolateName } from 'matterviz/colors'
import {
  bool_from_param,
  bool_url_entry,
  sync_url_params as sync_params,
  type UrlParamEntry,
} from 'svelte-widgets/url-params'

type PageState = Parameters<typeof replaceState>[1]
// color_scale param: valid d3 interpolate names, defaulting to Viridis
export const url_color_scale = {
  default: `interpolateViridis` as D3InterpolateName,
  read: (params: URLSearchParams): D3InterpolateName => {
    const value = params.get(`color_scale`) ?? ``
    return is_d3_interpolate_name(value) ? value : url_color_scale.default
  },
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
// Filter configuration shared by browser-history snapshots and table exports.
export type FilterConfig = {
  training: Record<string, TrainFilterMode>
  openness: readonly Openness[]
  targets?: Partial<Record<TargetOutput, TrainFilterMode>> // absent = default (require F)
  fs_mode?: FsMode
}
// minimal structural model shape keeps this module decoupled from $lib/types
type FilterableModel = {
  training_sets: string[]
  openness: Openness
  targets: string
}
const is_one_of = <Value extends string>(
  options: readonly Value[],
  value: unknown,
): value is Value => typeof value === `string` && options.includes(value as Value)

// Split a model.targets string like `EFS_GM` into its predicted outputs and the
// force/stress computation mode: prefix letters are E/F/S/H outputs, the suffix
// holds G(radient)/D(irect) plus M when the model also predicts magmoms.
// Exported so filter UIs can tally models per output with the same semantics.
export function parse_targets(targets: string): {
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
  show_selected_only = $state(false)
  private readonly training_entries = $derived(Object.entries(this.training))
  private readonly target_entries = $derived(Object.entries(this.targets))

  constructor(readonly training_sets: string[]) {}

  // number of active non-default constraints (drives filter-button badges)
  get n_active(): number {
    return (
      this.training_entries.length +
      (this.openness.length < OPENNESS_OPTIONS.length ? 1 : 0) +
      (this.targets_param === DEFAULT_TARGETS_PARAM ? 0 : 1)
    )
  }

  matches = (model: FilterableModel): boolean => this.matches_except(model, {})

  // Facet counts retain every other constraint, including other entries in the
  // same training/target filter, so they describe the current view.
  matches_except = (
    model: FilterableModel,
    ignored: {
      training?: string
      target?: TargetOutput
      openness?: boolean
      fs_mode?: boolean
    },
  ): boolean => {
    if (!ignored.openness && !this.openness.includes(model.openness)) return false
    const { outputs, fs_mode } = parse_targets(model.targets)
    const outputs_ok = this.target_entries.every(
      ([key, mode]) =>
        key === ignored.target || outputs.has(key) === (mode === `require`),
    )
    if (!outputs_ok) return false
    // direct/gradient also drops models without any force/stress prediction
    if (!ignored.fs_mode && this.fs_mode !== `any` && fs_mode !== this.fs_mode)
      return false
    return this.training_entries.every(
      ([key, mode]) =>
        key === ignored.training ||
        model.training_sets.includes(key) === (mode === `require`),
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

  apply = (config: FilterConfig): void => {
    // Snapshots survive deploys; discard obsolete constraints that no current
    // checkbox or URL parameter could represent.
    this.training = Object.fromEntries(
      Object.entries(config.training).filter(
        ([key, mode]) =>
          this.training_sets.includes(key) && is_one_of(TRAIN_FILTER_MODES, mode),
      ),
    )
    // filter OPENNESS_OPTIONS (not spread the config) to keep canonical order and
    // drop invalid tokens from stale snapshots
    const shown = OPENNESS_OPTIONS.filter((op) => config.openness.includes(op))
    this.openness = shown.length > 0 ? shown : [...OPENNESS_OPTIONS]
    this.targets = Object.fromEntries(
      Object.entries(config.targets ?? DEFAULT_TARGETS).filter(
        ([key, mode]) =>
          is_one_of(target_output_keys, key) && is_one_of(TRAIN_FILTER_MODES, mode),
      ),
    )
    this.fs_mode = is_one_of(FS_MODES, config.fs_mode) ? config.fs_mode : `any`
  }

  // Copy active filters so snapshots do not change with subsequent UI edits.
  get config(): FilterConfig {
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
    this.show_selected_only = bool_from_param(params, `selected_only`)
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
      bool_url_entry(`selected_only`, this.show_selected_only),
    ]
  }
}

export function sync_url_params(entries: UrlParamEntry[], state: PageState): void {
  sync_params(entries, location, (url) => replaceState(url, state))
}

let last_navigation: AfterNavigate | undefined
const popstate_queries = new WeakMap<Promise<void>, URLSearchParams>()

// Bind during component init. Wait for the router's first navigation before writing;
// late-mounted components read immediately, since Kit does not replay afterNavigate
// for them. Readers receive the latest navigation as well as the current query params.
export function bind_url_params(
  read_params: ((params: URLSearchParams, navigation: AfterNavigate) => void) | null,
  entries: () => UrlParamEntry[],
): void {
  let url_ready = $state(false)

  // Kit's history metadata retains the original router URL after shallow edits.
  // Capture the browser's restored query before destination components can write.
  onNavigate(({ type, complete }) => {
    if (type === `popstate` && !popstate_queries.has(complete)) {
      popstate_queries.set(complete, new URLSearchParams(location.search))
    }
  })
  const read_url = (
    navigation: AfterNavigate,
    params = popstate_queries.get(navigation.complete) ?? page.url.searchParams,
  ) => {
    last_navigation = navigation
    read_params?.(params, navigation)
    url_ready = true
  }
  afterNavigate(read_url)

  $effect(() => {
    const navigation = last_navigation
    // Normal navigations read the router's original query, before any URL writes.
    // Late mounts read shallow edits, which Kit's replaceState omits from page.url.
    if (!url_ready && navigation) {
      untrack(() => read_url(navigation, new URLSearchParams(location.search)))
    }
    if (!url_ready) return
    sync_url_params(entries(), page.state)
  })
}
