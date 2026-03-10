<script lang="ts">
  import { MetricScatter, MetricsTable, SelectToggle } from '$lib'
  import * as labels from '$lib/labels'
  import type { DiscoverySet } from '$lib/types'
  import HullConstructionNote from './hull-construction-note.md'

  let discovery_set: DiscoverySet = $state(`unique_prototypes`)
</script>

<h1>Crystal Stability Prediction Metrics</h1>

<SelectToggle
  bind:selected={discovery_set}
  options={Object.entries(labels.DISCOVERY_SET_LABELS).map(
    ([value, { label, description: tooltip, link }]) => ({
      value,
      label,
      tooltip,
      link,
    }),
  )}
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
  />
</section>

<HullConstructionNote />

<h2>F1 classification score vs model parameters</h2>

The F1 score is the harmonic mean of precision and recall. It is a measure of the model's
ability to correctly identify hypothetical crystals in the WBM test set as lying on or
below the Materials Project convex hull.

<MetricScatter
  x_prop={labels.HYPERPARAMS.model_params}
  y_prop={labels.ALL_METRICS.F1}
  style="height: 400px"
/>
