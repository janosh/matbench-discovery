<script lang="ts">
  import { MetricsTable, MODELS, SelectToggle } from '$lib'
  import { DynamicScatter } from '$lib/plot'
  import { scatter_axis_label } from '$lib/plot/DynamicScatter.svelte'
  import * as labels from '$lib/labels'
  import type { DiscoverySet } from '$lib/types'
  import HullConstructionNote from './hull-construction-note.md'

  let discovery_set: DiscoverySet = $state(`unique_prototypes`)
  let show_energy_only = $state(false)

  // axis selections for the model-comparison scatter, bound so the section title
  // tracks whatever properties the user picks
  let scatter_x = $state(labels.HYPERPARAMS.model_params.key)
  let scatter_y = $state(labels.ALL_METRICS.F1.key)
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
    show_energy_only_toggle
  />
</section>

<HullConstructionNote />

<h2>{@html scatter_axis_label(scatter_y)} vs {@html scatter_axis_label(scatter_x)}</h2>

The F1 score is the harmonic mean of precision and recall. It is a measure of the model's
ability to correctly identify hypothetical crystals in the WBM test set as lying on or
below the Materials Project convex hull. Use the axis/color/size selectors to compare
models across any pair of metrics and metadata.

<DynamicScatter
  models={MODELS}
  bind:x_key={scatter_x}
  bind:y_key={scatter_y}
  style="height: 800px"
/>
