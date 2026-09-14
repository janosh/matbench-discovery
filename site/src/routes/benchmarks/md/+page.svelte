<script lang="ts">
  import TestSet from '$lib/benchmark/TestSet.svelte'
  import TaskNavigation from '$lib/benchmark/TaskNavigation.svelte'
  import MetricsTable from '$lib/table/MetricsTable.svelte'
  import type { ModelData } from '$lib/types'
  import { ACTIVE_MODELS, make_table_filters } from '$lib/models.svelte'
  import {
    MD_METRICS,
    METADATA_COLS,
    scatter_axis_label,
    task_page_visible_cols,
  } from '$lib/labels'
  import { CMDS_CONFIG, DEFAULT_CMDS_CONFIG } from '$lib/combined-scores.svelte'
  import DynamicScatter from '$lib/plot/DynamicScatter.svelte'
  import RadarChart from '$lib/plot/RadarChart.svelte'
  import { bind_url_params } from '$lib/url-state.svelte'
  import MdNote from './md-note.md'

  // show only MD metrics and metadata columns
  const visible_cols = task_page_visible_cols(...Object.values(MD_METRICS))

  const has_md_metrics = (model: ModelData) => model.metrics?.md != null

  // Public cost-vs-fidelity view, sorted by CMDS (highest first).
  let plot = $state({
    x: MD_METRICS.md_run_time_sec.key,
    y: MD_METRICS.md_combined_score.key,
  })

  const filters = make_table_filters()

  bind_url_params(filters.read, () => filters.url_entries)
</script>

<h1 id="molecular-dynamics-metrics">
  Molecular Dynamics Metrics <span class="beta-badge">beta</span>
</h1>

<p>
  This task tests how well ML force fields reproduce structural, thermodynamic, and
  vibrational observables of ab-initio molecular dynamics (AIMD) at finite temperature.
</p>
<p>
  Reference data: <a href="#test-set">DynaMat</a>. The combined molecular dynamics score
  (CMDS, higher is better) combines trajectory fidelity with rollout speed.
</p>

<TaskNavigation />

<h2 id="leaderboard">Leaderboard</h2>
<p>
  <strong>Preliminary:</strong> the reference set and metrics are evolving, so rankings
  may change. Private-label energy and force errors are diagnostics, excluded from CMDS.
  {#if !ACTIVE_MODELS.some(has_md_metrics)}
    No models have reported MD metrics yet.
  {/if}
</p>

<section class="full-bleed">
  <MetricsTable
    model_filter={has_md_metrics}
    col_filter={(col) => visible_cols[col.key] ?? true}
    default_sort={{ column: MD_METRICS.md_combined_score.key, dir: `desc` }}
    {filters}
  />
</section>

<details style="margin-block: 1em">
  <summary>Adjust score weights</summary>
  <figure class="task-weights">
    <RadarChart
      size={260}
      config={CMDS_CONFIG}
      default_config={DEFAULT_CMDS_CONFIG}
      title_label={MD_METRICS.md_combined_score}
    />
    <figcaption>
      Drag the knob to reweight CMDS components; the table and plots update live. Hover
      the ⓘ icon for how CMDS is computed.
    </figcaption>
  </figure>
</details>

<h2 id="model-comparison" style="text-align: center">
  {@html scatter_axis_label(plot.y)} vs {@html scatter_axis_label(plot.x)}
</h2>
<p>
  This defaults to a cost-vs-fidelity Pareto: each model's total rollout wall time against
  its CMDS, with marker size showing model parameters and color the training-set size. Use
  the axis/color/size selectors to compare any pair of metrics: the RDF, ADF and vDOS
  errors range from 0% (perfect match with the AIMD reference) to 100% (as different from
  the reference as an ideal gas / non-overlapping distributions).
</p>

<DynamicScatter
  models={ACTIVE_MODELS}
  model_filter={has_md_metrics}
  bind:x_key={plot.x}
  bind:y_key={plot.y}
  color_key={METADATA_COLS.n_training_materials.key}
  show_pareto_frontier
  style="height: 800px"
/>

<TestSet task="md" />

<h2 id="methodology">Methodology</h2>
<details>
  <summary>Simulation protocol, reference data, and metric definitions</summary>
  <p>
    Each model runs NVT simulations from the same initial structures and thermodynamic
    conditions as the reference first-principles trajectories. The resulting trajectories
    are compared via radial distribution functions (RDF), angular distribution functions
    (ADF), pressure distributions from the stress tensor trace, and the vibrational
    density of states (vDOS) obtained from the velocity autocorrelation function.
    Energy-fluctuation and force RMSEs are shown as maintainer-computed private-label
    diagnostics when available, but they are excluded from CMDS.
  </p>
  <MdNote />
</details>

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
