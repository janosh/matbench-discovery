<script lang="ts">
  import { MetricScatter, MetricsTable, SelectToggle } from '$lib'
  import {
    DISCOVERY_METRICS,
    DISCOVERY_SET_LABELS,
    METADATA_COLS,
    METRICS,
  } from '$lib/labels'
  import type { DiscoverySet } from '$lib/types'

  let discovery_set: DiscoverySet = $state(`unique_prototypes`)
</script>

<h1>Crystal Stability Prediction Metrics</h1>

<SelectToggle
  bind:selected={discovery_set}
  options={Object.entries(DISCOVERY_SET_LABELS).map(
    ([value, { label, description: tooltip, link }]) => ({
      value,
      label,
      tooltip,
      link,
    }),
  )}
/>

<MetricsTable
  col_filter={(col) =>
    [...Object.values(DISCOVERY_METRICS), ...Object.values(METADATA_COLS)].includes(col)}
  {discovery_set}
  style="width: 100%;"
/>

<h3>F1 classification score over time</h3>

The F1 score is the harmonic mean of precision and recall. It is a measure of the model's
ability to correctly identify hypothetical crystals in the WBM test set as lying on or
below the Materials Project convex hull.

<MetricScatter
  x_prop={METADATA_COLS.date_added}
  y_prop={METRICS.F1}
  style="height: 400px;"
/>

<h3>F1 classification score vs model parameters</h3>

<MetricScatter
  x_prop={METADATA_COLS.model_params}
  y_prop={METRICS.F1}
  style="height: 400px;"
/>
