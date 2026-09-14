<script lang="ts">
  import TestSet from '$lib/benchmark/TestSet.svelte'
  import TaskNavigation from '$lib/benchmark/TaskNavigation.svelte'
  import MetricsTable from '$lib/table/MetricsTable.svelte'
  import DiscoverySetToggle from '$lib/DiscoverySetToggle.svelte'
  import { ACTIVE_MODELS, make_table_filters } from '$lib/models.svelte'
  import DynamicScatter from '$lib/plot/DynamicScatter.svelte'
  import { bind_url_params } from '$lib/url-state.svelte'
  import { valid_query_param } from 'svelte-widgets/url-params'
  import * as labels from '$lib/labels'
  import { DISCOVERY_SETS, type DiscoverySet } from '$lib/types'
  import HullConstructionNote from './hull-construction-note.md'

  const default_discovery_set: DiscoverySet = `unique_prototypes`
  // Axis selections also drive the model-comparison section title.
  let plot = $state({
    x: labels.HYPERPARAMS.model_params.key,
    y: labels.ALL_METRICS.F1.key,
  })
  const discovery_sets = new Set<DiscoverySet>(DISCOVERY_SETS)

  let discovery_set: DiscoverySet = $state(default_discovery_set)
  const filters = make_table_filters()
  let visible_models = $derived(
    ACTIVE_MODELS.filter(
      (model) =>
        model.metrics?.discovery?.[discovery_set] != null && filters.matches(model),
    ),
  )

  const read_url_params = (params: URLSearchParams) => {
    discovery_set = valid_query_param(
      params,
      `set`,
      default_discovery_set,
      discovery_sets,
    )
    filters.read(params)
  }
  bind_url_params(read_url_params, () => [
    [`set`, discovery_set, default_discovery_set],
    ...filters.url_entries,
  ])
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
    model_filter={(model) => model.metrics?.discovery?.[discovery_set] != null}
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
  models={visible_models}
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
