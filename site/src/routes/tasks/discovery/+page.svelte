<script lang="ts">
  import { MetricScatter, MetricsTable, MODELS, SelectToggle } from '$lib'
  import { DISCOVERY_METRICS, DISCOVERY_SET_LABELS, METADATA_COLS } from '$lib/metrics'
  import type { DiscoverySet } from '$lib/types'
  import type { Point } from 'elementari'

  // Default column visibility
  let visible_cols: Record<string, boolean> = $state({
    ...Object.fromEntries(
      [...DISCOVERY_METRICS, ...METADATA_COLS].map((col) => [col.label, true]),
    ),
    'κ<sub>SRME</sub>': false,
  })

  let discovery_set: DiscoverySet = $state(`unique_prototypes`)
  let f1_tooltip_point: Point | null = $state(null)
  let hovered = $state(false)

  let filtered_models = $derived(
    Object.values(MODELS).filter(
      (md) => md.metrics?.discovery?.[discovery_set]?.F1 != null,
    ),
  )
</script>

<h1>Crystal Stability Prediction Metrics</h1>

<SelectToggle
  bind:selected={discovery_set}
  options={Object.entries(DISCOVERY_SET_LABELS).map(
    ([value, { title, tooltip, link }]) => ({ value, label: title, tooltip, link }),
  )}
/>

<MetricsTable
  col_filter={(col) => visible_cols[col.label] ?? true}
  {discovery_set}
  style="width: 100%;"
/>

<h3>F1 classification score of models over time</h3>
<p>
  The F1 score is the harmonic mean of precision and recall. It is a measure of the
  model's ability to correctly identify hypothetical crystals in the WBM test set as lying
  on or below the Materials Project convex hull.
</p>

<MetricScatter
  models={filtered_models}
  metric="discovery.{discovery_set}.F1"
  y_label="F1 Score (higher better)"
  bind:tooltip_point={f1_tooltip_point}
  bind:hovered
  style="margin: 2em 0; width: 100%; height: 300px;"
/>

<style>
  h3 {
    text-align: center;
  }
</style>
