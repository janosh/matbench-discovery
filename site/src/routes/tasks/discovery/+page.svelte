<script lang="ts">
  import { MetricsTable, MODELS, SelectToggle } from '$lib'
  import { make_table_filters } from '$lib/models.svelte'
  import { DynamicScatter } from '$lib/plot'
  import type { SortState } from '$lib/url-state.svelte'
  import {
    bind_url_params,
    sort_from_query,
    sort_url_entries,
    valid_query_param,
  } from '$lib/url-state.svelte'
  import * as labels from '$lib/labels'
  import { DISCOVERY_SETS, type DiscoverySet, type ModelData } from '$lib/types'
  import HullConstructionNote from './hull-construction-note.md'

  const default_discovery_set: DiscoverySet = `unique_prototypes`
  const default_scatter_x = labels.HYPERPARAMS.model_params.key
  const default_scatter_y = labels.ALL_METRICS.F1.key
  const default_sort: SortState = { column: default_scatter_y, dir: `desc` }
  const discovery_sets = new Set<DiscoverySet>(DISCOVERY_SETS)

  let discovery_set: DiscoverySet = $state(default_discovery_set)
  const filters = make_table_filters()
  let sort = $state({ ...default_sort })
  let scatter_path_overrides = $derived(
    Object.fromEntries(
      Object.values(labels.DISCOVERY_METRICS).map(({ key }) => [
        key,
        `metrics.discovery.${discovery_set}.${key}`,
      ]),
    ),
  )
  const has_discovery_metrics = (model: ModelData): boolean => {
    const discovery = model.metrics?.discovery
    return typeof discovery === `object` && discovery?.[discovery_set] != null
  }
  let visible_models = $derived(
    MODELS.filter((model) => has_discovery_metrics(model) && filters.matches(model)),
  )

  // axis selections for the model-comparison scatter, bound so the section title
  // tracks whatever properties the user picks
  let scatter_x = $state(default_scatter_x)
  let scatter_y = $state(default_scatter_y)

  const scatter_keys = labels.scatter_options_by_key
  const read_url_params = (params: URLSearchParams) => {
    discovery_set = valid_query_param(
      params,
      `set`,
      default_discovery_set,
      discovery_sets,
    )
    filters.read(params)
    sort = sort_from_query(params, default_sort)
    scatter_x = valid_query_param(params, `x`, default_scatter_x, scatter_keys)
    scatter_y = valid_query_param(params, `y`, default_scatter_y, scatter_keys)
  }
  bind_url_params(read_url_params, () => [
    [`set`, discovery_set, default_discovery_set],
    ...filters.url_entries,
    ...sort_url_entries(sort, default_sort),
    [`x`, scatter_x, default_scatter_x],
    [`y`, scatter_y, default_scatter_y],
  ])
</script>

<h1>Crystal Stability Prediction Metrics</h1>

<div class="task-intro">
  <div>
    <p>
      This task measures how effectively a model can triage hypothetical WBM crystals for
      DFT validation. Models must identify structures that lie on or below a fixed
      DFT-computed Materials Project convex hull while minimizing costly false positives.
    </p>
    <p>
      Switch between the full test set, unique prototypes, and each model's 10,000
      most-stable predictions to compare overall accuracy, structural diversity, and
      fixed-budget discovery performance. See
      <a href="/tasks/discovery/tmi">Discovery TMI</a> for calibration, element-level, and error-distribution
      diagnostics.
    </p>
  </div>
</div>

<details class="methodology">
  <summary>Methodology: fixed DFT convex hull</summary>
  <HullConstructionNote />
</details>

<SelectToggle
  bind:selected={discovery_set}
  options={labels.discovery_set_toggle_options}
/>
<section class="full-bleed">
  <MetricsTable
    col_filter={(col) =>
      [
        labels.METADATA_COLS.model_name,
        ...Object.values(labels.DISCOVERY_METRICS),
        labels.METADATA_COLS.links,
        labels.METADATA_COLS.date_added,
      ].includes(col)}
    {discovery_set}
    model_filter={has_discovery_metrics}
    {filters}
    bind:sort
  />
</section>

<h2>
  {@html labels.scatter_axis_label(scatter_y)} vs {@html labels.scatter_axis_label(
    scatter_x,
  )}
</h2>

The F1 score is the harmonic mean of precision and recall. It is a measure of the model's
ability to correctly identify hypothetical crystals in the WBM test set as lying on or
below the Materials Project convex hull. Use the axis/color/size selectors to compare
models across any pair of metrics and metadata.

<!-- color by energy MAE: an orthogonal 3rd axis (the default color is F1, which is
already the y-axis here, so it wastes the color channel) -->
<DynamicScatter
  models={visible_models}
  bind:x_key={scatter_x}
  bind:y_key={scatter_y}
  color_key={labels.ALL_METRICS.MAE.key}
  label_path_overrides={scatter_path_overrides}
  style="height: 800px"
/>

<style>
  .methodology {
    margin-block: 1em;
  }
  .methodology summary {
    cursor: pointer;
    font-weight: 600;
  }
</style>
