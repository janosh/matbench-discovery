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
  } from '$lib/combined-scores.svelte'
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

  const has_md_metrics = (model: ModelData) =>
    model.metrics?.md != null && typeof model.metrics.md === `object`

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
  // default-sort by the combined MD score (CMDS), best (highest) first
  let sort = $state({ ...default_sort })

  const read_url_params = (params: URLSearchParams) => {
    scatter_x = valid_query_param(params, `x`, default_scatter_x, scatter_options_by_key)
    scatter_y = valid_query_param(params, `y`, default_scatter_y, scatter_options_by_key)
    sort = sort_from_query(params, default_sort)
    apply_weights_param(params.get(`weights`), CMDS_CONFIG, DEFAULT_CMDS_CONFIG)
  }
  bind_url_params(read_url_params, () => [
    [`x`, scatter_x, default_scatter_x],
    [`y`, scatter_y, default_scatter_y],
    [`sort`, sort.column, default_sort.column],
    [`dir`, sort.dir, default_sort.dir],
    // custom CMDS weights (vDOS,ADF,speed,pressure); omitted at defaults
    [`weights`, weights_to_param(CMDS_CONFIG, DEFAULT_CMDS_CONFIG)],
  ])
</script>

<h1>Molecular Dynamics Metrics <span class="beta-badge">beta</span></h1>

<MdNote />

<div class="intro bleed-1400">
  <p>
    This task evaluates how well machine-learning interatomic potentials reproduce
    structural, thermodynamic and vibrational observables of ab-initio molecular dynamics
    (AIMD) trajectories at finite temperature. Each model runs NVT simulations from the
    same initial structures and thermodynamic conditions as the reference first-principles
    trajectories. The resulting trajectories are compared via radial distribution
    functions (RDF), angular distribution functions (ADF), pressure distributions from the
    stress tensor trace, and the vibrational density of states (vDOS) obtained from the
    velocity autocorrelation function. Energy-fluctuation and force RMSEs are shown as
    maintainer-computed private-label diagnostics when available, but they are excluded
    from CMDS. This modeling task was introduced in
    <a href="https://arxiv.org/abs/2607.03433">arXiv:2607.03433</a>.
    {#if !MODELS.some(has_md_metrics)}
      No models have reported MD metrics yet.
    {/if}
  </p>
  <figure class="cmds-weights">
    <RadarChart
      size={260}
      config={CMDS_CONFIG}
      default_config={DEFAULT_CMDS_CONFIG}
      title_label={MD_METRICS.md_combined_score}
      on_change={(cfg) => update_models_cmds(MODELS, cfg as typeof CMDS_CONFIG)}
    />
    <figcaption>
      Drag the knob to reweight which CMDS components matter to you; the table and plots
      update live. Hover the ⓘ icon for how CMDS is computed.
    </figcaption>
  </figure>
</div>

<section class="full-bleed">
  <MetricsTable
    model_filter={has_md_metrics}
    col_filter={(col) => visible_cols[col.key] ?? true}
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
  .intro {
    display: flex;
    flex-wrap: wrap;
    gap: 1em 2em;
    align-items: center;
  }
  .intro > p {
    flex: 1 1 30em;
  }
  .cmds-weights {
    flex: 0 1 22em;
    margin: 0 auto;
  }
  .cmds-weights figcaption {
    margin-top: 1em;
    font-size: 0.85em;
    color: var(--text-muted, inherit);
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
