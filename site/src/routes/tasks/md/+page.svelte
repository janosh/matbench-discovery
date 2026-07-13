<script lang="ts">
  import { MetricsTable, type ModelData, MODELS } from '$lib'
  import {
    MD_METRICS,
    METADATA_COLS,
    scatter_axis_label,
    scatter_options_by_key,
    task_page_visible_cols,
  } from '$lib/labels'
  import type { SortState } from '$lib/url-state.svelte'
  import {
    CMDS_CONFIG,
    DEFAULT_CMDS_CONFIG,
    update_models_cmds,
  } from '$lib/combined-scores.svelte'
  import { DynamicScatter, RadarChart } from '$lib/plot'
  import { make_table_filters } from '$lib/models.svelte'
  import {
    apply_weights_param,
    bind_url_params,
    sort_from_query,
    sort_url_entries,
    valid_query_param,
    weights_to_param,
  } from '$lib/url-state.svelte'
  import MdNote from './md-note.md'

  // show only MD metrics and metadata columns
  const visible_cols = task_page_visible_cols(...Object.values(MD_METRICS))

  const has_md_metrics = (model: ModelData) => model.metrics?.md != null

  // headline MD view now that CMDS folds in rollout speed: a cost-vs-fidelity Pareto -
  // wall time (x) vs CMDS (y), with marker size = model params and color = training-set
  // size (the two scaling levers). All-public axes (force_rmse is a maintainer-only
  // diagnostic future public submissions lack, which would leave the plot mostly empty)
  const default_scatter_x = MD_METRICS.md_run_time_sec.key
  const default_scatter_y = MD_METRICS.md_combined_score.key
  const default_sort: SortState = {
    column: MD_METRICS.md_combined_score.key,
    dir: `desc`,
  }

  let scatter_x = $state(default_scatter_x)
  let scatter_y = $state(default_scatter_y)
  // default-sort by the combined MD score (CMDS), best (highest) first
  let sort = $state({ ...default_sort })
  const filters = make_table_filters()

  const read_url_params = (params: URLSearchParams) => {
    scatter_x = valid_query_param(params, `x`, default_scatter_x, scatter_options_by_key)
    scatter_y = valid_query_param(params, `y`, default_scatter_y, scatter_options_by_key)
    sort = sort_from_query(params, default_sort)
    filters.read(params)
    apply_weights_param(params.get(`weights`), CMDS_CONFIG, DEFAULT_CMDS_CONFIG)
  }
  bind_url_params(read_url_params, () => [
    [`x`, scatter_x, default_scatter_x],
    [`y`, scatter_y, default_scatter_y],
    ...sort_url_entries(sort, default_sort),
    ...filters.url_entries,
    // custom CMDS weights (vDOS,ADF,speed,pressure); omitted at defaults
    [`weights`, weights_to_param(CMDS_CONFIG, DEFAULT_CMDS_CONFIG)],
  ])
</script>

<h1>Molecular Dynamics Metrics <span class="beta-badge">beta</span></h1>

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
  {#if !MODELS.some(has_md_metrics)}
    No models have reported MD metrics yet.
  {/if}
</p>

<div class="task-intro bleed-1400">
  <!-- wrapper div: the markdown renders multiple top-level elements which would
  otherwise each become their own flex item -->
  <div><MdNote /></div>
  <figure class="task-weights">
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
    {filters}
  />
</section>

<h2>{@html scatter_axis_label(scatter_y)} vs {@html scatter_axis_label(scatter_x)}</h2>
<p>
  This defaults to a cost-vs-fidelity Pareto: each model's total rollout wall time against
  its CMDS, with marker size showing model parameters and color the training-set size. Use
  the axis/color/size selectors to compare any pair of metrics: the RDF, ADF and vDOS
  errors range from 0% (perfect match with the AIMD reference) to 100% (as different from
  the reference as an ideal gas / non-overlapping distributions).
</p>

<DynamicScatter
  models={MODELS}
  model_filter={has_md_metrics}
  bind:x_key={scatter_x}
  bind:y_key={scatter_y}
  color_key={METADATA_COLS.n_training_materials.key}
  show_pareto_frontier
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
