<script lang="ts">
  import kappa_103_analysis from '$figs/kappa-103-analysis.json.gz'
  import { by_date_added_desc, MetricsTable, type ModelData, MODELS } from '$lib'
  import { DynamicScatter, KappaParityPlot } from '$lib/plot'
  import { scatter_axis_label } from '$lib/plot/DynamicScatter.svelte'
  import { has_kappa_parity_model } from '$lib/kappa-parity'
  import { ALL_METRICS, HYPERPARAMS, METADATA_COLS } from '$lib/labels'
  import { get_nested_number } from '$lib/metrics'
  import { format_num } from 'matterviz'
  import KappaNote from './kappa-note.md'
  import { SvelteSet } from 'svelte/reactivity'
  import KappaSrmeScatter from './KappaSrmeScatter.svelte'
  import PhononFreqParity from './PhononFreqParity.svelte'
  import PhononRobustnessTable from './PhononRobustnessTable.svelte'

  // Default column visibility
  let visible_cols: Record<string, boolean> = $state({
    // Hide other metrics
    ...Object.fromEntries(
      Object.values(ALL_METRICS).map((col) => [col.label, false]),
    ),
    // Show all metadata
    ...Object.fromEntries(
      Object.values(METADATA_COLS).map((col) => [col.label, true]),
    ),
    // Show phonon metrics
    [ALL_METRICS.κ_SRME.label]: true,
    [ALL_METRICS.κ_SRE.label]: true,
  })

  // Models with generated kappa parity assets or per-material diagnostics, for the
  // DFT-vs-ML inspector below
  const diagnostics_keys = new SvelteSet(
    kappa_103_analysis.models.map((entry) => entry.key),
  )
  const kappa_models = MODELS.filter((model) =>
    has_kappa_parity_model(model.model_key) ||
    diagnostics_keys.has(model.model_key ?? ``)
  )
  let selected_key = $state(kappa_models[0]?.model_key)
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
  let sort_mode = $state<SortMode>(`kappa`)
  const srme_path = `${ALL_METRICS.κ_SRME.path}.${ALL_METRICS.κ_SRME.key}`
  const kappa_srme = (model: ModelData) => get_nested_number(model, srme_path) ?? Infinity
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
  // suffix each option with the value it's sorted by (κSRME or date; nothing for A–Z)
  const sort_label = (model: ModelData): string => {
    if (sort_mode === `date`) return ` (${model.date_added})`
    if (sort_mode === `name`) return ``
    const srme = get_nested_number(model, srme_path)
    return srme == null ? `` : ` (${format_num(srme, `.3~f`)})`
  }

  // axis selections for the model-comparison scatter, bound so the section title
  // tracks whatever properties the user picks
  let scatter_x = $state(HYPERPARAMS.model_params.key)
  let scatter_y = $state(ALL_METRICS.κ_SRME.key)
</script>

<h1>MLFF Phonon Modeling Metrics</h1>

<section class="full-bleed">
  <MetricsTable col_filter={(col) => visible_cols[col.label] ?? true} />
</section>

<KappaNote />

<h2>Failure Modes</h2>
<p>
  κ<sub>SRME</sub> assigns its maximum error of 2 to materials where the prediction
  pipeline breaks down entirely: imaginary phonon modes after ML relaxation (the model
  predicts an unstable structure), symmetry broken during relaxation, or a crashed κ
  calculation. This table shows how much of each model's κ<sub>SRME</sub> comes from
  such outright failures (and how many of those have imaginary modes as the known
  cause), alongside the Wasserstein-1 distance between ML and DFT phonon frequency
  spectra (a κ-independent measure of phonon accuracy that doesn't suffer from error
  compounding in the thermal conductivity calculation).
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
  >PBE reference data for the 103 MP structures</a> that form the test set for this task.
  Use the axis/color/size selectors to compare models across any pair of metrics and
  metadata. Clicking a point selects that model in the inspector below.
</p>

<DynamicScatter
  models={MODELS}
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
    <!-- {#key} rebuilds the native <select> on sort change; moving <option> nodes via a
    keyed {#each} alone doesn't repaint the native dropdown on some platforms -->
    {#key sort_mode}
      <select bind:value={selected_key}>
        {#each sorted_models as model (model.model_key)}
          <option value={model.model_key}>{model.model_name}{sort_label(model)}</option>
        {/each}
      </select>
    {/key}
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
      crystal system &mdash; failures concentrated at low/high κ or in specific
      symmetries point to systematic weaknesses. Hollow markers at κ<sub>SRME</sub> =
      2 are censored values (the κ calculation failed), not measurements. Right:
      quantile-quantile parity of the ML vs DFT phonon frequency spectra &mdash;
      points below the diagonal mean the model predicts too-soft phonons (under-stiff
      force constants), above means too-stiff.
    </p>
    <div class="diagnostics-grid bleed-1400">
      <KappaSrmeScatter
        entry={selected_diagnostics}
        base={kappa_103_analysis}
        style="height: 420px"
      />
      <PhononFreqParity entry={selected_diagnostics} style="height: 420px" />
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
  /* fixed width so the appended (value) in options doesn't resize the select */
  .kappa-model-select select:first-of-type {
    width: 18em;
    overflow: hidden;
    white-space: nowrap;
    text-overflow: ellipsis;
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
