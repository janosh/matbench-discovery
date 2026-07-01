<script lang="ts">
  import { afterNavigate } from '$app/navigation'
  import { page } from '$app/state'
  import { MetricsTable, type ModelData, MODELS } from '$lib'
  import {
    ALL_METRICS,
    MD_METRICS,
    METADATA_COLS,
    scatter_axis_label,
    scatter_options_by_key,
  } from '$lib/labels'
  import type { SortDir } from '$lib/types'
  import { DynamicScatter } from '$lib/plot'
  import { sort_from_query, sync_url_params, valid_query_param } from '$lib/url-state'
  import MdNote from './md-note.md'

  // show only MD metrics and metadata columns
  const visible_cols: Record<string, boolean> = Object.fromEntries([
    ...Object.entries(ALL_METRICS).map(([key, col]): [string, boolean] => [
      col.label,
      key in MD_METRICS,
    ]),
    ...Object.values(METADATA_COLS).map((col): [string, boolean] => [col.label, true]),
  ])

  // guard against null since typeof null === `object` (a null metrics.md is not data)
  const has_md_metrics = (model: ModelData) =>
    model.metrics?.md != null && typeof model.metrics.md === `object`
  const n_md_models = MODELS.filter(has_md_metrics).length

  const default_scatter_x = MD_METRICS.md_force_rmse.key
  const default_scatter_y = MD_METRICS.md_rdf_error.key
  const default_sort: { column: string; dir: SortDir } = {
    column: MD_METRICS.md_combined_score.key,
    dir: `desc`,
  }

  let scatter_x = $state(default_scatter_x)
  let scatter_y = $state(default_scatter_y)
  // default-sort by the combined MD score (CMDS), best (highest) first. matterviz sorts
  // by the column's key (falling back to label), so use the key, not the 'CMDS' label
  let sort = $state({ ...default_sort })
  let url_ready = $state(false)

  afterNavigate(() => {
    const params = page.url.searchParams
    scatter_x = valid_query_param(params, `x`, default_scatter_x, scatter_options_by_key)
    scatter_y = valid_query_param(params, `y`, default_scatter_y, scatter_options_by_key)
    sort = sort_from_query(params, default_sort)
    url_ready = true
  })

  $effect(() => {
    if (!url_ready) return
    sync_url_params(
      [
        [`x`, scatter_x, default_scatter_x],
        [`y`, scatter_y, default_scatter_y],
        [`sort`, sort.column, default_sort.column],
        [`dir`, sort.dir, default_sort.dir],
      ],
      page.state,
    )
  })
</script>

<h1>Molecular Dynamics Metrics <span class="beta-badge">beta</span></h1>

<MdNote />

<p>
  This task evaluates how well machine-learning interatomic potentials reproduce
  structural, thermodynamic and vibrational observables of ab-initio molecular dynamics
  (AIMD) trajectories at finite temperature. Each model runs NVT simulations from the same
  initial structures and thermodynamic conditions as the reference first-principles
  trajectories. The resulting trajectories are compared via radial distribution functions
  (RDF), angular distribution functions (ADF), pressure distributions from the stress
  tensor trace, and the vibrational density of states (vDOS) obtained from the velocity
  autocorrelation function. Single-point energy-fluctuation and force RMSEs on the
  reference frames complement these trajectory-level observables.
  {#if n_md_models === 0}
    No models have reported MD metrics yet.
  {/if}
</p>

<section class="full-bleed">
  <!-- ?? true: columns absent from visible_cols (e.g. the sticky model name) default to
  shown, matching the phonons and landing-page tables; visible_cols only hides non-MD metrics -->
  <MetricsTable
    model_filter={has_md_metrics}
    col_filter={(col) => visible_cols[col.label] ?? true}
    bind:sort
  />
</section>

<h2>{@html scatter_axis_label(scatter_y)} vs {@html scatter_axis_label(scatter_x)}</h2>
<p>
  RDF, ADF and vDOS errors range from 0% (perfect match with the AIMD reference) to 100%
  (as different from the reference as an ideal gas / non-overlapping distributions). Use
  the axis/color/size selectors to compare models across any pair of metrics and metadata.
</p>

<DynamicScatter
  models={MODELS}
  model_filter={has_md_metrics}
  bind:x_key={scatter_x}
  bind:y_key={scatter_y}
  color_key={MD_METRICS.md_combined_score.key}
  style="height: 800px"
/>

<style>
  .beta-badge {
    font-size: 0.45em;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.08em;
    vertical-align: middle;
    padding: 2px 7px;
    border-radius: 5px;
    color: orange;
    background: color-mix(in oklab, orange 18%, transparent);
    border: 1px solid color-mix(in oklab, orange 45%, transparent);
  }
</style>
