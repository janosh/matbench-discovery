<script lang="ts">
  import { by_date_added_desc, MetricsTable, type ModelData, MODELS } from '$lib'
  import { KappaParityPlot, MetricScatter } from '$lib/plot'
  import { has_kappa_parity_model } from '$lib/kappa-parity'
  import { ALL_METRICS, HYPERPARAMS, METADATA_COLS } from '$lib/labels'
  import { get_nested_number } from '$lib/metrics'
  import { format_num } from 'matterviz'
  import KappaNote from './kappa-note.md'

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

  // Models with generated kappa parity assets, for the DFT-vs-ML inspector below
  const kappa_models = MODELS.filter((model) => has_kappa_parity_model(model.model_key))
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
  // suffix each option with the value it's sorted by (κSRME or date; nothing for A–Z)
  const sort_label = (model: ModelData): string => {
    if (sort_mode === `date`) return ` (${model.date_added})`
    if (sort_mode === `name`) return ``
    const srme = get_nested_number(model, srme_path)
    return srme == null ? `` : ` (${format_num(srme, `.3~f`)})`
  }
</script>

<h1>MLFF Phonon Modeling Metrics</h1>

<section class="full-bleed">
  <MetricsTable col_filter={(col) => visible_cols[col.label] ?? true} />
</section>

<KappaNote />

<h2>κ<sub>SRME</sub> vs Model Parameters</h2>
<p>
  κ<sub>SRME</sub> ranges from 0 to 2, the lower the better. This metric was introduced in
  <a href="https://arxiv.org/abs/2408.00755v4">arXiv:2408.00755v4</a>. This modeling task
  would not have been possible without the
  <a href="https://github.com/atztogo/phonondb">PhononDB</a>
  and the help of Atsushi Togo who kindly shared the
  <a
    href="https://github.com/atztogo/phonondb/blob/main/README.md#url-links-to-phono3py-finite-displacement-method-inputs-of-103-compounds-on-mdr-at-nims-pbe"
  >PBE reference data for the 103 MP structures</a> that form the test set for this task.
</p>

<MetricScatter
  x_prop={HYPERPARAMS.model_params}
  y_prop={ALL_METRICS.κ_SRME}
  style="height: 400px"
  on_point_click={({ metadata }) => {
    const key = metadata?.model_key
    if (typeof key === `string` && has_kappa_parity_model(key)) selected_key = key
  }}
/>

{#if selected_model}
  <label class="kappa-model-select">
    Compare model:
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
</style>
