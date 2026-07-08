<script lang="ts">
  import kappa_103_analysis from '$figs/kappa-103-analysis.jsonl'
  import { by_date_added_desc, MetricsTable, ModelSelect, MODELS } from '$lib'
  import type { ModelData } from '$lib'
  import { DynamicScatter, KappaParityPlot } from '$lib/plot'
  import { make_table_filters } from '$lib/models.svelte'
  import { DEFAULT_TABLE_SORT } from '$lib/table/MetricsTable.svelte'
  import {
    bind_url_params,
    sort_from_query,
    sort_url_entries,
    valid_query_param,
  } from '$lib/url-state.svelte'
  import { has_kappa_parity_model } from '$lib/parity/kappa-parity'
  import {
    ALL_METRICS,
    HYPERPARAMS,
    scatter_axis_label,
    scatter_options_by_key,
    task_page_visible_cols,
  } from '$lib/labels'
  import { get_nested_number } from '$lib/metrics'
  import { format_num } from 'matterviz'
  import KappaNote from './kappa-note.md'
  import { SvelteSet } from 'svelte/reactivity'
  import KappaSrmeScatter from './KappaSrmeScatter.svelte'
  import PhononFreqParity from './PhononFreqParity.svelte'
  import PhononRobustnessTable from './PhononRobustnessTable.svelte'

  // Default column visibility: metadata + phonon metrics only
  const visible_cols = task_page_visible_cols(ALL_METRICS.κ_SRME, ALL_METRICS.κ_SRE)

  // Models with generated kappa parity assets or per-material diagnostics, for the
  // DFT-vs-ML inspector below
  const diagnostics_keys = new SvelteSet(
    kappa_103_analysis.models.map((entry) => entry.key),
  )
  const kappa_srme_path = `${ALL_METRICS.κ_SRME.path}.${ALL_METRICS.κ_SRME.key}`
  const kappa_sre_path = `${ALL_METRICS.κ_SRE.path}.${ALL_METRICS.κ_SRE.key}`
  const has_phonon_metrics = (model: ModelData): boolean =>
    get_nested_number(model, kappa_srme_path) != null &&
    get_nested_number(model, kappa_sre_path) != null
  const kappa_models = MODELS.filter(
    (model) =>
      has_kappa_parity_model(model.model_key) ||
      diagnostics_keys.has(model.model_key ?? ``),
  )
  const default_selected_key = kappa_models[0]?.model_key ?? ``
  let selected_key = $state(default_selected_key)
  let selected_model = $derived(
    kappa_models.find((model) => model.model_key === selected_key),
  )

  // sort the Compare-model dropdown by κ_SRME (best first), name, or submission date
  type SortMode = `kappa` | `name` | `date`
  const sort_options: { mode: SortMode; label: string }[] = [
    { mode: `kappa`, label: `κSRME` },
    { mode: `name`, label: `A–Z` },
    { mode: `date`, label: `date added` },
  ]
  const default_sort_mode: SortMode = `kappa`
  const default_scatter_x = HYPERPARAMS.model_params.key
  const default_scatter_y = ALL_METRICS.κ_SRME.key
  const sort_modes = new Set(sort_options.map(({ mode }) => mode))
  const model_keys = new Set(kappa_models.flatMap((model) => model.model_key ?? []))

  let sort_mode = $state<SortMode>(default_sort_mode)
  const kappa_srme = (model: ModelData) =>
    get_nested_number(model, kappa_srme_path) ?? Infinity
  const sort_compare: Record<SortMode, (m1: ModelData, m2: ModelData) => number> = {
    kappa: (m1, m2) => kappa_srme(m1) - kappa_srme(m2),
    name: (m1, m2) => m1.model_name.localeCompare(m2.model_name),
    date: by_date_added_desc,
  }
  let sorted_models = $derived(kappa_models.toSorted(sort_compare[sort_mode]))
  // per-material diagnostics (SRME scatter + frequency parity) for the selected model
  let selected_diagnostics = $derived(
    kappa_103_analysis.models.find((entry) => entry.key === selected_key),
  )
  // model dropdown options, each suffixed with the value it's sorted by (κSRME or date;
  // nothing for A–Z); value carries the model_key so the bound option maps cleanly onto
  // selected_key (which URL state and scatter clicks also set)
  let model_options = $derived(
    sorted_models.map((model) => {
      const srme = get_nested_number(model, kappa_srme_path)
      const suffix = {
        name: ``,
        date: ` (${model.date_added})`,
        kappa: srme == null ? `` : ` (${format_num(srme, `.3~f`)})`,
      }[sort_mode]
      return { label: `${model.model_name}${suffix}`, value: model.model_key ?? `` }
    }),
  )

  // axis selections for the model-comparison scatter, bound so the section title
  // tracks whatever properties the user picks
  let scatter_x = $state(default_scatter_x)
  let scatter_y = $state(default_scatter_y)

  let table_sort = $state({ ...DEFAULT_TABLE_SORT })
  const filters = make_table_filters()

  const read_url_params = (params: URLSearchParams) => {
    selected_key = valid_query_param(params, `model`, default_selected_key, model_keys)
    sort_mode = valid_query_param(params, `model_sort`, default_sort_mode, sort_modes)
    scatter_x = valid_query_param(params, `x`, default_scatter_x, scatter_options_by_key)
    scatter_y = valid_query_param(params, `y`, default_scatter_y, scatter_options_by_key)
    table_sort = sort_from_query(params, DEFAULT_TABLE_SORT)
    filters.read(params)
  }
  bind_url_params(read_url_params, () => [
    [`model`, selected_key, default_selected_key],
    [`model_sort`, sort_mode, default_sort_mode],
    [`x`, scatter_x, default_scatter_x],
    [`y`, scatter_y, default_scatter_y],
    ...sort_url_entries(table_sort, DEFAULT_TABLE_SORT),
    ...filters.url_entries,
  ])
</script>

<h1>MLFF Phonon Modeling Metrics</h1>

<section class="full-bleed">
  <MetricsTable
    model_filter={has_phonon_metrics}
    col_filter={(col) => visible_cols[col.key] ?? true}
    bind:sort={table_sort}
    {filters}
  />
</section>

<KappaNote />

<h2>Failure Modes</h2>
<p>
  κ<sub>SRME</sub> assigns its maximum error of 2 to materials where the prediction
  pipeline breaks down entirely: imaginary phonon modes after ML relaxation (the model
  predicts an unstable structure), symmetry broken during relaxation, or a crashed κ
  calculation. This table shows how much of each model's κ<sub>SRME</sub> comes from such outright
  failures (and how many of those have imaginary modes as the known cause), alongside the Wasserstein-1
  distance between ML and DFT phonon frequency spectra (a κ-independent measure of phonon accuracy
  that doesn't suffer from error compounding in the thermal conductivity calculation).
</p>
<section class="full-bleed robustness-table">
  <PhononRobustnessTable />
</section>

<h2>{@html scatter_axis_label(scatter_y)} vs {@html scatter_axis_label(scatter_x)}</h2>
<p>
  κ<sub>SRME</sub> ranges from 0 to 2, the lower the better. This metric was introduced in
  <a href="https://arxiv.org/abs/2408.00755v4">arXiv:2408.00755v4</a>. This modeling task
  would not have been possible without the
  <a href="https://github.com/atztogo/phonondb">PhononDB</a>
  and the help of Atsushi Togo who kindly shared the
  <a
    href="https://github.com/atztogo/phonondb/blob/main/README.md#url-links-to-phono3py-finite-displacement-method-inputs-of-103-compounds-on-mdr-at-nims-pbe"
    >PBE reference data for the 103 MP structures</a
  > that form the test set for this task. Use the axis/color/size selectors to compare models
  across any pair of metrics and metadata. Clicking a point selects that model in the inspector
  below.
</p>

<DynamicScatter
  models={MODELS}
  model_filter={has_phonon_metrics}
  bind:x_key={scatter_x}
  bind:y_key={scatter_y}
  style="height: 800px"
  point_events={{
    onclick: ({ point }) => {
      const key = point.metadata?.model_key
      if (typeof key === `string` && kappa_models.some((mdl) => mdl.model_key === key)) {
        selected_key = key
      }
    },
  }}
/>

{#if selected_model}
  <label class="kappa-model-select">
    View model:
    <ModelSelect
      options={model_options}
      minSelect={1}
      maxSelect={1}
      style="width: 22em; border: 1px solid var(--border)"
      bind:value={
        () => model_options.find((opt) => opt.value === selected_key) ?? null,
        // Array.isArray only narrows the type; with maxSelect=1 value is never an array
        (option) => {
          if (option && !Array.isArray(option)) selected_key = option.value
        }
      }
    />
    <select bind:value={sort_mode} aria-label="Sort models by">
      {#each sort_options as opt (opt.mode)}
        <option value={opt.mode}>sort: {opt.label}</option>
      {/each}
    </select>
  </label>
  <KappaParityPlot model={selected_model} />

  {#if selected_diagnostics}
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
  .kappa-model-select {
    display: flex;
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
    grid-template-columns: 1fr 1fr;
    gap: 1em;
  }
  @media (max-width: 900px) {
    .diagnostics-grid {
      grid-template-columns: 1fr;
    }
  }
</style>
