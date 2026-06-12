<script lang="ts">
  import { MetricsTable, type ModelData, MODELS } from '$lib'
  import { ALL_METRICS, MD_METRICS, METADATA_COLS } from '$lib/labels'
  import { DynamicScatter } from '$lib/plot'
  import { scatter_axis_label } from '$lib/plot/DynamicScatter.svelte'

  // show only MD metrics and metadata columns
  const visible_cols: Record<string, boolean> = {
    ...Object.fromEntries(
      Object.values(ALL_METRICS).map((col) => [col.label, false]),
    ),
    ...Object.fromEntries(
      Object.values(METADATA_COLS).map((col) => [col.label, true]),
    ),
    ...Object.fromEntries(
      Object.values(MD_METRICS).map((col) => [col.label, true]),
    ),
  }

  const has_md_metrics = (model: ModelData) => typeof model.metrics?.md === `object`
  const n_md_models = MODELS.filter(has_md_metrics).length

  let scatter_x = $state(MD_METRICS.MD_force_RMSE.key)
  let scatter_y = $state(MD_METRICS.MD_RDF_error.key)
</script>

<h1>Molecular Dynamics Metrics</h1>

<p>
  This task evaluates how well machine-learning interatomic potentials reproduce
  structural, thermodynamic and vibrational observables of ab-initio molecular
  dynamics (AIMD) trajectories at finite temperature. Each model runs NVT simulations
  from the same initial structures and thermodynamic conditions as the reference
  first-principles trajectories. The resulting trajectories are compared via radial
  distribution functions (RDF), pressure distributions from the stress tensor trace,
  and the vibrational density of states (VDOS) obtained from the velocity
  autocorrelation function. Single-point energy and force RMSEs on the reference
  frames complement these trajectory-level observables.
  {#if n_md_models === 0}
    No models have reported MD metrics yet.
  {/if}
</p>

<section class="full-bleed">
  <MetricsTable col_filter={(col) => visible_cols[col.label] ?? true} />
</section>

<h2>{@html scatter_axis_label(scatter_y)} vs {@html scatter_axis_label(scatter_x)}</h2>
<p>
  RDF and VDOS errors range from 0% (perfect match with the AIMD reference) to 100%
  (as different from the reference as an ideal gas / non-overlapping spectra). Use the
  axis/color/size selectors to compare models across any pair of metrics and metadata.
</p>

<DynamicScatter
  models={MODELS}
  model_filter={has_md_metrics}
  bind:x_key={scatter_x}
  bind:y_key={scatter_y}
  color_key={MD_METRICS.MD_combined_error.key}
  style="height: 800px"
/>
