import { page } from '$app/state'
import { DATASETS } from '$lib'
import { ALL_METRICS, HYPERPARAMS, METADATA_COLS } from '$lib/labels'
import {
  discovery_task_tooltips,
  get_nested_value,
  is_finite_num,
  label_data_path,
  openness_tooltips,
  targets_tooltips,
} from '$lib/metrics'
import { ACTIVE_MODELS, MODELS } from '$lib/models.svelte'
import type { Author, Label, ModelData } from '$lib/types'
import { bind_url_params } from '$lib/url-state.svelte'
import { format_num } from 'matterviz'
import type { RowData } from 'matterviz/table'
import { SvelteSet } from 'svelte/reactivity'

const model_by_key = new Map(MODELS.map((model) => [model.model_key, model]))

// Site-wide set of models picked for side-by-side comparison (any number); insertion
// order is display order. Lives at module level like CPS_CONFIG so tables, cards, model
// pages and the comparison dialog all share it without prop drilling.
class Comparison {
  keys = new SvelteSet<string>()
  open = $state(false)

  get models(): ModelData[] {
    return [...this.keys].flatMap((key) => model_by_key.get(key) ?? [])
  }
  toggle = (key: string): void => {
    if (!this.keys.delete(key)) this.keys.add(key)
  }
  // "Compare with…" entry points: make sure the model is selected, then open the dialog
  open_with = (key: string): void => {
    this.keys.add(key)
    this.open = true
  }
  // replace the selection (unknown keys dropped), skipping the write when nothing changed
  set = (keys: string[]): void => {
    const next = keys.filter((key) => model_by_key.has(key))
    if (next.join(',') === [...this.keys].join(',')) return
    this.keys.clear()
    for (const key of next) this.keys.add(key)
  }
}
export const comparison = new Comparison()

// Leaderboard table hooks shared by MetricsTable and GeoOptMetricsTable: double-clicking
// a row toggles its model; `show_only` drops non-compared rows, else they get highlighted.
export const toggle_row_model = (event: Event, row: RowData): void => {
  if (typeof row.model_key !== `string`) return
  event.preventDefault()
  comparison.toggle(row.model_key)
}
export function mark_compared_rows<R extends { model_key: string; class?: string }>(
  rows: R[],
  show_only: boolean,
): R[] {
  return rows
    .filter((row) => !show_only || comparison.keys.has(row.model_key))
    .map((row) => {
      row.class =
        !show_only && comparison.keys.has(row.model_key) ? `highlight` : undefined
      return row
    })
}

// Mirror the selection into a `compare` URL param (comma-joined model keys) so comparisons
// are shareable and survive reloads. Unlike other URL-bound state, an absent param keeps
// the in-memory selection (it follows the user across pages) and is re-applied to every
// navigated-to URL. Call once during layout init.
export function bind_comparison_url(): void {
  bind_url_params(
    (params, { type }) => {
      const param = params.get(`compare`)
      if (param !== null) comparison.set(param.split(`,`))
      // a shared link opens the dialog; in-app navigation (e.g. a model link inside the
      // dialog) closes it
      comparison.open = type === `enter` && comparison.keys.size > 1
    },
    () => {
      void page.url.pathname // re-sync after navigations whose target URL lacked the param
      return [[`compare`, [...comparison.keys].join(`,`)]]
    },
  )
}

// A piece of a text cell: plain text, optionally linked (site-relative or external) and/or
// with a hover tooltip carrying details that don't fit in the cell
export type CellPart = { text: string; href?: string; title?: string }
// A comparison row is a metric/metadata label, optionally with a custom accessor for
// values not addressable by a dotted path (derived costs, protocol summaries). `parts`
// renders text cells richly (links, tooltips); `title` adds a tooltip to numeric cells.
export type CompareRow = Label & {
  value?: (model: ModelData) => unknown
  parts?: (model: ModelData) => CellPart[] | undefined
  title?: (model: ModelData) => string | undefined
}
export type CompareGroup = { title: string; href?: string; rows: CompareRow[] }

export const row_value = (row: CompareRow, model: ModelData): unknown =>
  row.value ? row.value(model) : get_nested_value(model, label_data_path(row))

const metrics = (...keys: (keyof typeof ALL_METRICS)[]): CompareRow[] =>
  keys.map((key) => ALL_METRICS[key])
// text rows sort/filter on their concatenated parts; `parts` returning nothing hides the cell
const text_row = (
  key: string,
  label: string,
  description: string,
  parts: (model: ModelData) => CellPart[] | undefined,
): CompareRow => ({
  key,
  label,
  description,
  parts,
  value: (model) =>
    parts(model)
      ?.map((part) => part.text)
      .join(``),
})
// interleave parts with a plain separator, e.g. `MPtrj + OMat24`
const join_parts = (parts: CellPart[], separator: string): CellPart[] =>
  parts.flatMap((part, idx) => (idx > 0 ? [{ text: separator }, part] : [part]))
// person links in order of preference; affiliation as hover detail
const person_part = ({ name, url, github, orcid, affiliation }: Author): CellPart => ({
  text: name,
  href: url ?? github ?? orcid,
  title: affiliation,
})
const MAX_AUTHORS = 3

// device-hours per hardware entry plus how the total was estimated
const training_cost_title = ({ training_cost }: ModelData): string | undefined =>
  training_cost &&
  [
    ...training_cost.entries.map(
      ({ count, hardware, hours_per_device }) =>
        `${count}× ${hardware} × ${format_num(hours_per_device, `.3~s`)} h`,
    ),
    training_cost.reason,
  ]
    .filter(Boolean)
    .join(`<br>`)

// Cost rows double as the x-axis menu of the cost-vs-accuracy scatter in
// ModelComparison.svelte, hence lower=better on all of them.
export const COST_ROWS: CompareRow[] = [
  {
    ...HYPERPARAMS.model_params,
    better: `lower`,
    title: ({ n_estimators = 1, model_params }) =>
      n_estimators > 1
        ? `Ensemble of ${n_estimators} models with ${format_num(model_params, `.3~s`)} parameters each`
        : undefined,
  },
  { ...METADATA_COLS.n_training_materials, better: `lower` },
  { ...METADATA_COLS.n_training_structures, better: `lower` },
  {
    key: `training_gpu_hours`,
    label: `Training compute`,
    unit: `GPU·h`,
    format: `.3~s`,
    better: `lower`,
    description: `Reported training cost summed over all listed devices (device count × hours per device)`,
    value: (model) =>
      model.training_cost?.entries.reduce(
        (sum, { count, hours_per_device }) => sum + count * hours_per_device,
        0,
      ),
    title: training_cost_title,
  },
  { ...ALL_METRICS.md_run_time_sec, label: `MD speed` },
  { ...ALL_METRICS.md_max_gpu_mem_gb, label: `MD VRAM` },
  { ...ALL_METRICS.diatomics_run_time_sec, label: `Diatomics speed` },
]

// Rows grouped by task; rows where no compared model has a value are hidden in the UI.
export const COMPARE_GROUPS: CompareGroup[] = [
  {
    title: `Model`,
    rows: [
      ALL_METRICS.CPS,
      text_row(
        `targets`,
        `Targets`,
        `Quantities the model predicts and how forces/stress are obtained`,
        ({ targets }) => [{ text: targets, title: targets_tooltips[targets] }],
      ),
      text_row(
        `openness`,
        `Openness`,
        `Whether source code and training data are open`,
        ({ openness }) => [{ text: openness, title: openness_tooltips[openness] }],
      ),
      text_row(
        `training_sets`,
        `Training data`,
        `Datasets the model was trained on (linked to their data pages)`,
        ({ training_sets }) =>
          join_parts(
            training_sets.map((key) => {
              const { name, slug, n_structures, n_materials } = DATASETS[key]
              const from = n_materials ? ` from ${format_num(n_materials)} materials` : ``
              return {
                text: key,
                href: `/data/${slug}`,
                title: `${name}: ${format_num(n_structures)} structures${from}`,
              }
            }),
            ` + `,
          ),
      ),
      text_row(
        `authors`,
        `Authors`,
        `Model authors, linked to their website, GitHub or ORCID where available`,
        ({ authors }) => {
          const shown = authors.slice(0, MAX_AUTHORS).map(person_part)
          const rest = authors.slice(MAX_AUTHORS)
          if (rest.length > 0) {
            const text = `+${rest.length} more`
            shown.push({ text, title: rest.map((person) => person.name).join(`, `) })
          }
          return join_parts(shown, `, `)
        },
      ),
      text_row(
        `links`,
        `Links`,
        `Paper, code repository, documentation, checkpoint and package resources`,
        (model) =>
          join_parts(
            (
              [
                [`Paper`, model.paper],
                [`DOI`, model.doi],
                [`Repo`, model.repo],
                [`Docs`, model.docs],
                [`Checkpoint`, model.checkpoint_url],
                [`PyPI`, model.pypi],
              ] as const
            ).flatMap(([text, href]) => (href ? [{ text, href, title: href }] : [])),
            ` · `,
          ),
      ),
      text_row(
        `license`,
        `License`,
        `Licenses of the model code and checkpoint`,
        ({ license }) =>
          join_parts(
            [
              {
                text: license.code,
                href: license.code_url ?? undefined,
                title: `Code license${license.code_url_reason ? `: ${license.code_url_reason}` : ``}`,
              },
              {
                text: license.checkpoint,
                href: license.checkpoint_url ?? undefined,
                title: `Checkpoint license${license.checkpoint_url_reason ? `: ${license.checkpoint_url_reason}` : ``}`,
              },
            ],
            ` / `,
          ),
      ),
      {
        ...METADATA_COLS.benchmark_added,
        description: `Date the model was included on the benchmark leaderboard (linked to the pull request that added it)`,
        parts: ({ dates: { benchmark_added }, pr_url }) =>
          benchmark_added
            ? [
                {
                  text: benchmark_added,
                  href: pr_url ?? undefined,
                  title: pr_url
                    ? `Open the pull request that added this model`
                    : undefined,
                },
              ]
            : [],
      },
    ],
  },
  { title: `Cost`, rows: COST_ROWS },
  {
    title: `Discovery`,
    href: `/tasks/discovery`,
    rows: [
      ...metrics(`F1`, `DAF`, `Precision`, `Accuracy`, `MAE`, `RMSE`, `R2`),
      text_row(
        `test_task`,
        `Test task`,
        `How WBM stability was predicted: IS2RE-SR relaxes structures with the model first, IS2RE/IS2E predict from initial structures`,
        ({ test_task }) => [
          { text: test_task, title: discovery_task_tooltips[test_task] },
        ],
      ),
    ],
  },
  {
    title: `Geometry optimization`,
    href: `/tasks/geo-opt`,
    rows: [
      ...metrics(`RMSD`, `symmetry_match_1e-2`, `symmetry_decrease_1e-2`),
      text_row(
        `relaxation`,
        `Relaxation`,
        `ASE optimizer, force threshold and step limit used to relax structures`,
        (model) => {
          const { ase_optimizer, max_force, max_steps } =
            model.hyperparams?.evaluation ?? {}
          if (!ase_optimizer) return []
          const text = `${ase_optimizer} · f_max ${max_force ?? `?`} eV/Å · ${max_steps ?? `?`} steps`
          return [{ text }]
        },
      ),
    ],
  },
  {
    title: `Phonons`,
    href: `/tasks/phonons`,
    rows: metrics(`κ_SRME`, `κ_SRE`, `κ_failure_rate`, `imaginary_mode_rate`),
  },
  {
    title: `Molecular dynamics`,
    href: `/tasks/md`,
    rows: [
      ...metrics(
        `md_combined_score`,
        `md_vdos_error`,
        `md_adf_error`,
        `md_pressure_error`,
      ),
      text_row(
        `md_hardware`,
        `Hardware`,
        `Device the MD rollouts were timed on`,
        (model) => {
          const hardware = model.metrics?.md?.hardware
          return hardware ? [{ text: hardware }] : []
        },
      ),
    ],
  },
  {
    title: `Diatomics`,
    href: `/tasks/diatomics`,
    rows: [
      ...metrics(
        `diatomics_combined_score`,
        `pbe_energy_mae`,
        `pbe_force_mae`,
        `energy_jump`,
        `force_flips`,
      ),
      text_row(
        `diatomics_hardware`,
        `Hardware`,
        `Device the diatomic sweep was timed on`,
        (model) => {
          const hardware = model.metrics?.diatomics?.hardware
          return hardware ? [{ text: hardware }] : []
        },
      ),
    ],
  },
]

export type CompareCell = {
  text: string
  parts?: CellPart[] // rich rendering of a text cell (links, per-part tooltips)
  title?: string // hover tooltip for a numeric cell
  rank?: number // 1-based competition rank among `field` models with a value
  n?: number // number of field models with a value
  best?: boolean // best among the compared models (ties share it)
}

// One cell per compared model: formatted value plus leaderboard rank in the active cohort
// and whether it's the best of the compared set. Text rows get neither.
export function compare_cells(
  row: CompareRow,
  models: readonly ModelData[],
  field: readonly ModelData[] = ACTIVE_MODELS,
): CompareCell[] {
  const values = models.map((model) => row_value(row, model))
  const sign = row.better === `lower` ? -1 : 1 // larger sign * value is always better
  const field_nums = row.better
    ? field.map((model) => row_value(row, model)).filter(is_finite_num)
    : []
  const numbers = values.filter(is_finite_num)
  const best = numbers.length > 1 ? Math.max(...numbers.map((num) => sign * num)) : NaN
  return values.map((value, idx) => {
    const model = models[idx]
    if (!is_finite_num(value)) {
      const parts = row.parts?.(model) ?? []
      const text = parts.length > 0 ? parts.map((part) => part.text).join(``) : value
      if (typeof text !== `string` || !text) return { text: `–` }
      return parts.length > 0 ? { text, parts } : { text }
    }
    const text = format_num(value, row.format ?? `.3~g`)
    const title = row.title?.(model)
    const cell: CompareCell = title ? { text, title } : { text }
    if (!row.better) return cell
    const rank = field_nums.filter((other) => sign * other > sign * value).length + 1
    return { ...cell, rank, n: field_nums.length, best: sign * value === best }
  })
}
