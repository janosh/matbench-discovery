<script lang="ts">
  import { MetricsTable, type ModelData, MODELS } from '$lib'
  import {
    MD_METRICS,
    scatter_axis_label,
    scatter_options_by_key,
    task_page_visible_cols,
  } from '$lib/labels'
  import type { SortDir } from '$lib/types'
  import {
    CMDS_CONFIG,
    DEFAULT_CMDS_CONFIG,
    update_models_cmds,
  } from '$lib/md_combined_score.svelte'
  import { DynamicScatter, RadarChart } from '$lib/plot'
  import {
    apply_weights_param,
    bind_url_params,
    sort_from_query,
    valid_query_param,
    weights_to_param,
  } from '$lib/url-state.svelte'
  import MdNote from './md-note.md'

  // show only MD metrics and metadata columns
  const visible_cols = task_page_visible_cols(...Object.values(MD_METRICS))

  // guard against null since typeof null === `object` (a null metrics.md is not data)
  const has_md_metrics = (model: ModelData) =>
    model.metrics?.md != null && typeof model.metrics.md === `object`
  const n_md_models = MODELS.filter(has_md_metrics).length

  // default to public observables: force_rmse is a maintainer-only diagnostic that
  // future public submissions won't have, which would leave the scatter mostly empty
  const default_scatter_x = MD_METRICS.md_pressure_error.key
  const default_scatter_y = MD_METRICS.md_vdos_error.key
  const default_sort: { column: string; dir: SortDir } = {
    column: MD_METRICS.md_combined_score.key,
    dir: `desc`,
  }

  let scatter_x = $state(default_scatter_x)
  let scatter_y = $state(default_scatter_y)
  // default-sort by the combined MD score (CMDS), best (highest) first. matterviz sorts
  // by the column's key (falling back to label), so use the key, not the 'CMDS' label
  let sort = $state({ ...default_sort })

  const read_url_params = (params: URLSearchParams) => {
    scatter_x = valid_query_param(params, `x`, default_scatter_x, scatter_options_by_key)
    scatter_y = valid_query_param(params, `y`, default_scatter_y, scatter_options_by_key)
    sort = sort_from_query(params, default_sort)
    apply_weights_param(params.get(`weights`), CMDS_CONFIG)
  }
  bind_url_params(read_url_params, () => [
    [`x`, scatter_x, default_scatter_x],
    [`y`, scatter_y, default_scatter_y],
    [`sort`, sort.column, default_sort.column],
    [`dir`, sort.dir, default_sort.dir],
    // custom CMDS weights (ADF,vDOS,pressure); omitted at defaults
    [`weights`, weights_to_param(CMDS_CONFIG, DEFAULT_CMDS_CONFIG)],
  ])
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
  autocorrelation function. Energy-fluctuation and force RMSEs are shown as
  maintainer-computed private-label diagnostics when available, but they are excluded from
  CMDS.
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

<figure class="cmds-weights">
  <figcaption>
    CMDS is computed on the fly as a weighted mean of the component subscores (1 −
    error/100), not stored with model submissions. Drag the knob to reweight which
    observables matter to you; the table and plots update live.
  </figcaption>
  <RadarChart
    size={260}
    config={CMDS_CONFIG}
    default_config={DEFAULT_CMDS_CONFIG}
    title_label={MD_METRICS.md_combined_score}
    on_change={(cfg) => update_models_cmds(MODELS, cfg as typeof CMDS_CONFIG)}
  />
</figure>

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
  .cmds-weights {
    display: flex;
    flex-wrap: wrap;
    gap: 1.5em;
    align-items: center;
    margin: 1em auto;
    max-width: 45em;
  }
  .cmds-weights figcaption {
    font-size: 0.92em;
    color: var(--text-color-muted, inherit);
  }
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
