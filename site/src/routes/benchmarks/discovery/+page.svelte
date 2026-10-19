<script lang="ts">
  import TestSet from '#lib/benchmark/TestSet.svelte'
  import TaskNavigation from '#lib/benchmark/TaskNavigation.svelte'
  import MetricsTable from '#lib/table/MetricsTable.svelte'
  import DiscoverySetToggle from '#lib/DiscoverySetToggle.svelte'
  import { ACTIVE_MODELS, task_page_filters } from '#lib/models.svelte.js'
  import DynamicScatter from '#lib/plot/DynamicScatter.svelte'
  import { valid_query_param } from 'svelte-widgets/url-params'
  import * as labels from '#lib/labels.js'
  import { DISCOVERY_SETS, type DiscoverySet, type ModelData } from '#lib/types.js'
  import HullConstructionNote from './hull-construction-note.md'

  const default_discovery_set: DiscoverySet = `unique_prototypes`
  // Axis selections also drive the model-comparison section title.
  let plot = $state({
    x: labels.HYPERPARAMS.model_params.key,
    y: labels.ALL_METRICS.F1.key,
  })
  const discovery_sets = new Set<DiscoverySet>(DISCOVERY_SETS)

  let discovery_set: DiscoverySet = $state(default_discovery_set)
  const has_discovery_metrics = (model: ModelData) =>
    model.metrics?.discovery?.[discovery_set] != null
  const { filters, in_view } = task_page_filters(
    has_discovery_metrics,
    (params) => {
      discovery_set = valid_query_param(
        params,
        `set`,
        default_discovery_set,
        discovery_sets,
      )
    },
    () => [[`set`, discovery_set, default_discovery_set]],
  )
</script>

<h1 id="crystal-stability-prediction-metrics">Crystal Stability Prediction Metrics</h1>

<p>
  This task measures how effectively a model can triage hypothetical crystals for DFT
  validation. Models must identify structures that lie on or below a fixed DFT-computed
  Materials Project convex hull while minimizing costly false positives.
</p>
<p>
  Reference data: the <a href="#test-set">WBM test set</a>. F1 (higher is better) is the
  harmonic mean of precision and recall. Switch between the full test set and unique
  prototypes to compare overall accuracy against performance on the deduplicated, more
  out-of-distribution subset.
</p>

<TaskNavigation />

<h2 id="leaderboard">Leaderboard</h2>
<DiscoverySetToggle bind:selected={discovery_set} />
<section class="full-bleed">
  <MetricsTable
    col_filter={(col) =>
      [
        labels.METADATA_COLS.model_name,
        ...Object.values(labels.DISCOVERY_METRICS),
        labels.METADATA_COLS.links,
        labels.METADATA_COLS.benchmark_added,
      ].includes(col)}
    {discovery_set}
    model_filter={has_discovery_metrics}
    {filters}
    default_sort={{ column: labels.ALL_METRICS.F1.key, dir: `desc` }}
  />
</section>

<h2 id="model-comparison" style="text-align: center">
  {@html labels.scatter_axis_label(plot.y)} vs {@html labels.scatter_axis_label(plot.x)}
</h2>

<p>
  Use the axis/color/size selectors to compare models across any pair of metrics and
  metadata. The default view plots F1 against model size, colored by energy MAE.
</p>

<!-- color by energy MAE: an orthogonal 3rd axis (the default color is F1, which is
already the y-axis here, so it wastes the color channel) -->
<DynamicScatter
  models={ACTIVE_MODELS}
  model_filter={in_view}
  bind:x_key={plot.x}
  bind:y_key={plot.y}
  color_key={labels.ALL_METRICS.MAE.key}
  {discovery_set}
  style="height: 800px"
/>

<h2 id="diagnostics">Diagnostics</h2>
<p>
  Explore <a href="/benchmarks/discovery/tmi"
    >calibration, element-level errors, and error distributions</a
  > to see where stability predictions succeed or fail.
</p>

<TestSet task="discovery" />

<h2 id="methodology">Methodology</h2>
<HullConstructionNote />
