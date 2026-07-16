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
    PHONON_METRICS,
    scatter_axis_label,
    scatter_options_by_key,
    task_page_visible_cols,
  } from '$lib/labels'
  import { format_num } from 'matterviz'
  import KappaSrmeScatter from './KappaSrmeScatter.svelte'
  import PhononFreqParity from './PhononFreqParity.svelte'

  // Default column visibility: metadata + phonon metrics only
  const visible_cols = task_page_visible_cols(...Object.values(PHONON_METRICS))

  const kappa_srme = (model: ModelData): number =>
    model.metrics?.phonons?.kappa_103?.κ_SRME ?? Infinity
  const has_phonon_metrics = (model: ModelData): boolean =>
    Number.isFinite(kappa_srme(model))
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
  const default_scatter_x = PHONON_METRICS.κ_SRE.key
  const default_scatter_y = PHONON_METRICS.κ_SRME.key
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
      const srme = model.metrics?.phonons?.kappa_103?.κ_SRME
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
  const github_src_url = `https://github.com/janosh/matbench-discovery/blob/main/matbench_discovery`
</script>

<h1>MLFF Phonon Modeling Metrics</h1>

<div class="task-intro">
  <div>
    <p>
      This benchmark evaluates 300 K lattice thermal conductivity for 103 Materials
      Project crystals from the <a href={phonon_url}>PhononDB PBE test set</a>. After one
      simultaneous cell-and-site relaxation, finite-displacement forces yield second- and
      third-order force constants and each material's conductivity.
    </p>
    <p>
      κ<sub>SRME</sub> compares mode contributions before summation, while κ<sub>SRE</sub>
      compares only the final scalar conductivity, where cancellation can hide errors. Both
      range from 0 (perfect) to 2 (maximum error), with lower values better. κ<sub
        >SRD</sub
      > retains the sign of the scalar error to show systematic under- or overprediction.
    </p>
  </div>
</div>

<blockquote class="task-note">
  κ<sub>SRME</sub> follows the method of
  <a href="https://arxiv.org/abs/2408.00755v4">Póta et al.</a>. See the implementations
  for
  <a href="{github_src_url}/phonons/thermal_conductivity.py">
    phonon and conductivity prediction</a
  >
  and
  <a href="{github_src_url}/metrics/phonons.py">metric evaluation</a>.
</blockquote>

<h2>Leaderboard</h2>
<p>
  κ failed is the fraction of predictions whose κ<sub>SRME</sub> was censored to 2 because
  of imaginary modes, broken symmetry, or invalid conductivity data. A valid κ<sub
    >SRME</sub
  > of 2 is not a failure; Im(ω) separately reports the imaginary-mode rate.
</p>
<section class="full-bleed">
  <MetricsTable
    model_filter={has_phonon_metrics}
    col_filter={(col) => visible_cols[col.key] ?? true}
    bind:sort={table_sort}
    {filters}
  />
</section>

<h2>
  Model Comparison: {@html scatter_axis_label(scatter_y)} vs {@html scatter_axis_label(
    scatter_x,
  )}
</h2>
<p>
  This defaults to κ<sub>SRME</sub> versus κ<sub>SRE</sub>: stronger models sit toward the
  lower left, while separation between the axes reveals cancellation across mode
  contributions. Use the selectors to compare other properties; clicking a point opens
  that model below.
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
      Left: per-material κ<sub>SRME</sub> against DFT conductivity, colored by crystal
      system. Hollow markers are censored failures; valid κ<sub>SRME</sub> = 2 points stay filled.
      Right: quantile-quantile parity of the ML vs DFT phonon frequency spectra, colored by
      W₁(ω); points below the diagonal are too soft and points above are too stiff.
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
