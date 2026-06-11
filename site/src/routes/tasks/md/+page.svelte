<script lang="ts">
  import { MetricsTable, MODELS, type ModelData } from '$lib'
  import { DynamicScatter } from '$lib/plot'
  import { scatter_axis_label } from '$lib/plot/DynamicScatter.svelte'
  import { ALL_METRICS, MD_METRICS, METADATA_COLS } from '$lib/labels'
  import { get_nested_value } from '$lib/metrics'

  const visible_cols = [
    METADATA_COLS.model_name,
    ...Object.values(MD_METRICS),
    METADATA_COLS.links,
    METADATA_COLS.date_added,
  ]

  const has_md_metrics = (model: ModelData) =>
    typeof get_nested_value(model, `metrics.md`) === `object`

  let scatter_x = $state(ALL_METRICS.MD_force_RMSE.key)
  let scatter_y = $state(ALL_METRICS.MD_RDF_error.key)
</script>

<h1>Molecular Dynamics Metrics</h1>

<section class="full-bleed">
  <MetricsTable
    col_filter={(col) => visible_cols.includes(col)}
    show_non_compliant
  />
</section>

<h2>{@html scatter_axis_label(scatter_y)} vs {@html scatter_axis_label(scatter_x)}</h2>

<DynamicScatter
  models={MODELS}
  model_filter={has_md_metrics}
  bind:x_key={scatter_x}
  bind:y_key={scatter_y}
  color_key={ALL_METRICS.MD_combined_error.key}
  style="height: 800px"
/>
