<script lang="ts">
  import type { ModelData, ModelMetadata } from '$lib'
  import { ScatterPlot, type Point } from 'elementari'

  interface Props {
    models: ModelData[]
    metric: string
    y_label: string
    tooltip_point?: Point | null
    hovered?: boolean
    style?: string
    date_range?: [Date | null, Date | null]
    model_filter?: (model: ModelData) => boolean
    y_lim?: [number, number] | null
    color?: string | null
    [key: string]: unknown
  }

  let {
    models,
    metric,
    y_label,
    tooltip_point = $bindable(null),
    hovered = $bindable(false),
    style = ``,
    date_range = [null, null],
    model_filter = () => true,
    y_lim = null,
    color = null,
    ...rest
  }: Props = $props()

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
  let filtered_models = $derived(
    models
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
      }),
  )

  let y = $derived(
    filtered_models.map((model) => get_nested_val(model.metrics, metric) ?? 0),
  )
  let x = $derived(
    filtered_models.map((model) => new Date(model.date_added ?? 0).getTime()),
  )

  // Find active model by matching either the index or timestamp
  let active_model = $derived(
    tooltip_point
      ? filtered_models.find((model, idx) => {
          const model_time = new Date(model.date_added ?? 0).getTime()
          // Check if x matches either the index or timestamp
          return (
            idx + 1 === tooltip_point?.x || // index-based match
            Math.abs(model_time - (tooltip_point?.x ?? 0)) < 24 * 60 * 60 * 1000 // date-based match (within 1 day)
          )
        })
      : null,
  )

  // default fill color based on metric type
  const default_color_map: Record<string, string> = {
    'phonons.kappa_103.κ_SRME': `#ff6b6b`,
    'discovery.unique_prototypes.F1': `#ffbb54`,
    combined_performance_score: `#a78bfa`,
  }
  const fill_color = $derived(color ?? default_color_map[metric] ?? `#4dabf7`)

  // Create the data series with proper typing
  let point_style = $derived({
    fill: fill_color,
    radius: 5,
    stroke: `white`,
    stroke_width: 0.5,
  })
</script>

<ScatterPlot
  series={[{ x, y, point_style }]}
  x_format="%b %y"
  {y_label}
  bind:tooltip_point
  bind:hovered
  {style}
  y_lim={y_lim ?? undefined}
  {...rest}
>
  {#snippet tooltip({ x_formatted, y_formatted })}
    <div style="min-width: 10em;">
      {#if active_model}
        <strong>{active_model.model_name}</strong>
        <br />
        {@html y_label} = {@html y_formatted}
        <br />
        Added: {x_formatted}
      {/if}
    </div>
  {/snippet}
</ScatterPlot>
