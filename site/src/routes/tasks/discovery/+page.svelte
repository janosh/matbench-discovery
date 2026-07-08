<script lang="ts">
  import { MetricsTable, MODELS, SelectToggle } from '$lib'
  import { make_table_filters } from '$lib/models.svelte'
  import { DynamicScatter } from '$lib/plot'
  import { DEFAULT_TABLE_SORT } from '$lib/table/MetricsTable.svelte'
  import {
    bind_url_params,
    bool_from_param,
    bool_url_entry,
    sort_from_query,
    sort_url_entries,
    valid_query_param,
  } from '$lib/url-state.svelte'
  import * as labels from '$lib/labels'
  import { DISCOVERY_SETS, type DiscoverySet } from '$lib/types'
  import HullConstructionNote from './hull-construction-note.md'

  const default_discovery_set: DiscoverySet = `unique_prototypes`
  const default_scatter_x = labels.HYPERPARAMS.model_params.key
  const default_scatter_y = labels.ALL_METRICS.F1.key
  const discovery_sets = new Set<DiscoverySet>(DISCOVERY_SETS)

  let discovery_set: DiscoverySet = $state(default_discovery_set)
  let show_energy_only = $state(false)
  const filters = make_table_filters()
  let sort = $state({ ...DEFAULT_TABLE_SORT })

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
    show_energy_only = bool_from_param(params, `energy_only`)
    filters.read(params)
    sort = sort_from_query(params, DEFAULT_TABLE_SORT)
    scatter_x = valid_query_param(params, `x`, default_scatter_x, scatter_keys)
    scatter_y = valid_query_param(params, `y`, default_scatter_y, scatter_keys)
  }
  bind_url_params(read_url_params, () => [
    [`set`, discovery_set, default_discovery_set],
    bool_url_entry(`energy_only`, show_energy_only),
    ...filters.url_entries,
    ...sort_url_entries(sort, DEFAULT_TABLE_SORT),
    [`x`, scatter_x, default_scatter_x],
    [`y`, scatter_y, default_scatter_y],
  ])
</script>

<h1>Crystal Stability Prediction Metrics</h1>

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
    bind:show_energy_only
    {filters}
    bind:sort
    show_energy_only_toggle
  />
</section>

<HullConstructionNote />

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
  models={MODELS}
  bind:x_key={scatter_x}
  bind:y_key={scatter_y}
  color_key={labels.ALL_METRICS.MAE.key}
  style="height: 800px"
/>
