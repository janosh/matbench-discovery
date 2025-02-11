<script lang="ts">
  import type { ModelData } from '$lib'
  import { ScatterPlot } from 'elementari'
  import type { ModelMetadata } from './model-schema'

  export let models: ModelData[]
  export let metric: string
  export let y_label: string
  export let tooltip_point: { x: number; y: number } | null = null
  export let hovered: boolean = false
  export let style: string = ``
  export let date_range: [Date | null, Date | null] = [null, null]
  export let model_filter: (model: ModelData) => boolean = () => true

  // Add date range state
  const now = new Date()
  const ms_per_day = 24 * 60 * 60 * 1000
  const n_days_ago = new Date(now.getTime() - 180 * ms_per_day)

  // access nested values using a dotted path string
  function get_nested_val(obj: ModelMetadata, path: string) {
    return path.split(`.`).reduce((acc, key) => acc?.[key], obj)
  }

  // Filter out models where the metric is undefined and sort by date
  // Also filter by date range
  $: filtered_models = models
    .filter((model) => {
      const model_date = new Date(model.date_added ?? 0)
      return (
        get_nested_val(model.metrics, metric) &&
        model_date >= (date_range[0] ?? n_days_ago) &&
        model_date <= (date_range[1] ?? now) &&
        model_filter(model)
      )
    })
    .sort((model1, model2) => {
      const date_a = new Date(model1.date_added ?? 0)
      const date_b = new Date(model2.date_added ?? 0)
      return date_a.getTime() - date_b.getTime()
    })

  $: y = filtered_models.map((model) => get_nested_val(model.metrics, metric) ?? 0)
  $: x = filtered_models.map((model) => new Date(model.date_added ?? 0).getTime())

  // Find active model by matching either the index or timestamp
  $: active_model = tooltip_point
    ? filtered_models.find((model, idx) => {
        const model_time = new Date(model.date_added ?? 0).getTime()
        // Check if x matches either the index or timestamp
        return (
          idx + 1 === tooltip_point?.x || // index-based match
          Math.abs(model_time - (tooltip_point?.x ?? 0)) < 24 * 60 * 60 * 1000 // date-based match (within 1 day)
        )
      })
    : null
</script>

<ScatterPlot
  {x}
  {y}
  x_format="%b %y"
  {y_label}
  bind:tooltip_point
  bind:hovered
  {style}
  y_lim={[0, 2]}
  {...$$restProps}
>
  <div slot="tooltip" let:x_formatted let:y_formatted style="min-width: 10em;">
    {#if active_model}
      <strong>{active_model.model_name}</strong>
      <br />
      {y_label} = {y_formatted}
      <br />
      Added: {x_formatted}
    {/if}
  </div>
</ScatterPlot>
