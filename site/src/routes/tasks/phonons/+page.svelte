<script lang="ts">
  import kappa_103_analysis from '$figs/kappa-103-analysis.jsonl'
  import { by_benchmark_added_desc, MetricsTable, ModelSelect, ACTIVE_MODELS } from '$lib'
  import type { ModelData } from '$lib'
  import { DynamicScatter, KappaParityPlot } from '$lib/plot'
  import { make_table_filters } from '$lib/models.svelte'
  import type { SortState } from '$lib/url-state.svelte'
  import {
    bind_url_params,
    sort_from_query,
    sort_url_entries,
    valid_query_param,
  } from '$lib/url-state.svelte'
  import {
    ALL_METRICS,
    scatter_axis_label,
    scatter_options_by_key,
    task_page_visible_cols,
  } from '$lib/labels'
  import { get_nested_number, label_data_path } from '$lib/metrics'
  import { format_num } from 'matterviz'
  import KappaNote from './kappa-note.md'
  import KappaSrmeScatter from './KappaSrmeScatter.svelte'
  import PhononFreqParity from './PhononFreqParity.svelte'
  import PhononRobustnessTable from './PhononRobustnessTable.svelte'

  // Default column visibility: metadata + phonon metrics only
  const visible_cols = task_page_visible_cols(ALL_METRICS.κ_SRME, ALL_METRICS.κ_SRE)

  const kappa_srme_path = label_data_path(ALL_METRICS.κ_SRME)
  const kappa_sre_path = label_data_path(ALL_METRICS.κ_SRE)
  const has_phonon_metrics = (model: ModelData): boolean =>
    get_nested_number(model, kappa_srme_path) != null &&
    get_nested_number(model, kappa_sre_path) != null
  const kappa_srme = (model: ModelData): number =>
    get_nested_number(model, kappa_srme_path) ?? Infinity
  const leaderboard_models = ACTIVE_MODELS.filter(has_phonon_metrics)
  const default_selected_key =
    leaderboard_models.toSorted(
      (model_1, model_2) => kappa_srme(model_1) - kappa_srme(model_2),
    )[0]?.model_key ?? ``
  let selected_key = $state(default_selected_key)
  let selected_model = $derived(
    leaderboard_models.find((model) => model.model_key === selected_key),
  )

  // sort the Compare-model dropdown by κ_SRME (best first), name, or benchmark date
  type SortMode = `kappa` | `name` | `date`
  const sort_options: { mode: SortMode; label: string }[] = [
    { mode: `kappa`, label: `κSRME` },
    { mode: `name`, label: `A–Z` },
    { mode: `date`, label: `date added` },
  ]
  const default_sort_mode: SortMode = `kappa`
  const default_scatter_x = ALL_METRICS.κ_SRE.key
  const default_scatter_y = ALL_METRICS.κ_SRME.key
  const default_table_sort: SortState = {
    column: default_scatter_y,
    dir: `asc`,
  }
  const sort_modes = new Set(sort_options.map(({ mode }) => mode))
  const model_keys = new Set(leaderboard_models.map((model) => model.model_key))

  let sort_mode = $state<SortMode>(default_sort_mode)
  const sort_compare: Record<
    SortMode,
    (model_1: ModelData, model_2: ModelData) => number
  > = {
    kappa: (model_1, model_2) => kappa_srme(model_1) - kappa_srme(model_2),
    name: (model_1, model_2) => model_1.model_name.localeCompare(model_2.model_name),
    date: by_benchmark_added_desc,
  }
  let sorted_models = $derived(leaderboard_models.toSorted(sort_compare[sort_mode]))
  // per-material diagnostics (SRME scatter + frequency parity) for the selected model
  let selected_diagnostics = $derived(
    kappa_103_analysis.models.find((entry) => entry.key === selected_key),
  )
  // Model options include the value currently used for sorting (except A–Z).
  let model_options = $derived(
    sorted_models.map((model) => {
      const srme = get_nested_number(model, kappa_srme_path)
      const suffix = {
        name: ``,
        date: ` (${model.dates.benchmark_added})`,
        kappa: srme == null ? `` : ` (${format_num(srme, `.3~f`)})`,
      }[sort_mode]
      return {
        label: `${model.model_name}${suffix}`,
        value: model.model_key,
      }
    }),
  )

  // axis selections for the model-comparison scatter, bound so the section title
  // tracks whatever properties the user picks
  let scatter_x = $state(default_scatter_x)
  let scatter_y = $state(default_scatter_y)

  let table_sort = $state({ ...default_table_sort })
  const filters = make_table_filters()

  const read_url_params = (params: URLSearchParams) => {
    selected_key = valid_query_param(params, `model`, default_selected_key, model_keys)
    sort_mode = valid_query_param(params, `model_sort`, default_sort_mode, sort_modes)
    scatter_x = valid_query_param(params, `x`, default_scatter_x, scatter_options_by_key)
    scatter_y = valid_query_param(params, `y`, default_scatter_y, scatter_options_by_key)
    table_sort = sort_from_query(params, default_table_sort)
    filters.read(params)
  }
  bind_url_params(read_url_params, () => [
    [`model`, selected_key, default_selected_key],
    [`model_sort`, sort_mode, default_sort_mode],
    [`x`, scatter_x, default_scatter_x],
    [`y`, scatter_y, default_scatter_y],
    ...sort_url_entries(table_sort, default_table_sort),
    ...filters.url_entries,
  ])

  const phonon_url = `https://github.com/atztogo/phonondb/blob/main/README.md#url-links-to-phono3py-finite-displacement-method-inputs-of-103-compounds-on-mdr-at-nims-pbe`
</script>

<h1>MLFF Phonon Modeling Metrics</h1>

<div class="task-intro">
  <div>
    <p>
      This task evaluates whether machine-learning force fields reproduce lattice thermal
      conductivity at 300 K for 103 Materials Project crystals in the
      <a href="https://github.com/atztogo/phonondb">PhononDB</a> PBE test set. After one simultaneous
      cell-and-site relaxation, finite-displacement forces yield second- and third-order force
      constants, phonons, and each material's scalar conductivity.
    </p>
    <p>
      κ<sub>SRME</sub> (symmetric relative mean error) compares the individual phonon-mode
      contributions before they are summed, whereas κ<sub>SRE</sub>
      (symmetric relative error) compares only the final scalar conductivity. Both are on [0,
      2]: 0 is perfect, 2 is the maximum error, and <strong>lower is better</strong>.
      Over- and underpredicted mode contributions can cancel in the total, so a model can
      have a small κ<sub>SRE</sub> while retaining a larger κ<sub>SRME</sub>—the right
      total conductivity for the wrong microscopic reasons.
    </p>
  </div>
</div>

<!-- wrapper div: the markdown renders multiple top-level blockquotes -->
<div class="task-note"><KappaNote /></div>

<h2>Leaderboard</h2>
<p>
  The leaderboard requires both aggregate κ<sub>SRME</sub> and κ<sub>SRE</sub> values in model
  metadata; the inspector below provides the corresponding per-material views.
</p>
<section class="full-bleed">
  <MetricsTable
    model_filter={has_phonon_metrics}
    col_filter={(col) => visible_cols[col.key] ?? true}
    bind:sort={table_sort}
    {filters}
  />
</section>

<h2>Failure Diagnostics</h2>
<p>
  κ<sub>SRME</sub> assigns its maximum error of 2 to materials where the prediction
  pipeline breaks down entirely: imaginary phonon modes after ML relaxation (the model
  predicts an unstable structure), symmetry broken during relaxation, or a crashed κ
  calculation. This table shows how much of each model's κ<sub>SRME</sub> comes from such outright
  failures (and how many of those have imaginary modes as the known cause), alongside the Wasserstein-1
  distance between ML and DFT phonon frequency spectra (a κ-independent measure of phonon accuracy
  that doesn't suffer from error compounding in the thermal conductivity calculation). This
  section covers models with per-material analysis assets.
</p>
<section class="full-bleed robustness-table">
  <PhononRobustnessTable />
</section>

<h2>
  Model Comparison: {@html scatter_axis_label(scatter_y)} vs {@html scatter_axis_label(
    scatter_x,
  )}
</h2>
<p>
  This defaults to mode-resolved error against scalar-conductivity error for the same
  two-metric cohort as the leaderboard. Both default axes are lower-is-better, so stronger
  models sit toward the lower left; separation between them exposes cancellation across
  mode contributions. Use the axis/color/size selectors to compare other metrics and
  metadata. Clicking a point selects that model in the inspector below. The κ<sub
    >SRME</sub
  >
  metric was introduced in
  <a href="https://arxiv.org/abs/2408.00755v4">arXiv:2408.00755v4</a>; Atsushi Togo kindly
  shared the
  <a href={phonon_url}>PBE reference data for the 103 test structures</a>.
</p>

<DynamicScatter
  models={ACTIVE_MODELS}
  model_filter={has_phonon_metrics}
  bind:x_key={scatter_x}
  bind:y_key={scatter_y}
  style="height: 800px"
  point_events={{
    onclick: ({ point }) => {
      const model_key = point.metadata?.model_key
      if (typeof model_key === `string` && model_keys.has(model_key)) {
        selected_key = model_key
      }
    },
  }}
/>

<h2>Model Inspector</h2>
<p>Use the picker to inspect conductivity parity and per-material phonon errors.</p>
{#if selected_model}
  <label class="kappa-model-select">
    View model:
    <ModelSelect
      options={model_options}
      minSelect={1}
      maxSelect={1}
      style="width: 32em; max-width: 100%; border: 1px solid var(--border)"
      bind:selected={
        () => model_options.filter((option) => option.value === selected_key),
        (selected_options) => {
          const selected_option = selected_options[0]
          if (selected_option) selected_key = String(selected_option.value)
        }
      }
    />
    <select bind:value={sort_mode} aria-label="Sort models by">
      {#each sort_options as option (option.mode)}
        <option value={option.mode}>sort: {option.label}</option>
      {/each}
    </select>
  </label>

  {#if selected_diagnostics}
    <KappaParityPlot model={selected_model} />
    <h3 class="diagnostics-heading">
      Per-material κ<sub>SRME</sub> and phonon spectrum parity
    </h3>
    <p>
      Left: each material's κ<sub>SRME</sub> against its DFT conductivity, colored by
      crystal system &mdash; failures concentrated at low/high κ or in specific symmetries
      point to systematic weaknesses. Hollow markers at κ<sub>SRME</sub> = 2 are censored values
      (the κ calculation failed), not measurements. Right: quantile-quantile parity of the ML
      vs DFT phonon frequency spectra, colored by each material's spectrum W1 error. Tooltips
      include crystal system and quantile; points below the diagonal mean too-soft phonons,
      above means too-stiff.
    </p>
    <div class="diagnostics-grid bleed-1400">
      <KappaSrmeScatter
        entry={selected_diagnostics}
        base={kappa_103_analysis}
        style="height: 420px"
      />
      <PhononFreqParity
        entry={selected_diagnostics}
        base={kappa_103_analysis}
        style="height: 420px"
      />
    </div>
  {/if}
{/if}

<style>
  .task-note {
    margin-block: 1em 2em;
  }
  .kappa-model-select {
    display: flex;
    flex-wrap: wrap;
    gap: 0.5em;
    align-items: center;
    justify-content: center;
    margin-top: 1em;
  }
  .diagnostics-heading {
    text-align: center;
    margin-top: 1.5em;
  }
  .diagnostics-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(min(100%, 25rem), 1fr));
    gap: 1em;
  }
</style>
