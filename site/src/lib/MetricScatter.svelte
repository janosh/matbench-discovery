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

  // access nested values using a dotted path string
  function get_nested_val(obj: ModelMetadata, path: string) {
    return path.split(`.`).reduce((acc, key) => acc?.[key], obj)
  }

  // Filter out models where the metric is undefined and sort by date
  $: filtered_models = models
    .filter((model) => ![null, undefined].includes(get_nested_val(model.metrics, metric)))
    .sort((a, b) => {
      const date_a = new Date(a.date_added ?? 0)
      const date_b = new Date(b.date_added ?? 0)
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
  {...$$props}
  bind:tooltip_point
  bind:hovered
  {style}
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
