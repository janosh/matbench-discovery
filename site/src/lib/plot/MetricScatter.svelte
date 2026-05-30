<script lang="ts">
  import type { ModelData } from '$lib'
  import { MODELS } from '$lib'
  import { get_nested_value } from '$lib/metrics'
  import type { Label } from '$lib/types'
  import { type DataSeries, type PointStyle, ScatterPlot } from 'matterviz'
  import type { ComponentProps } from 'svelte'
  import type { HTMLAttributes } from 'svelte/elements'

  let {
    x_prop,
    y_prop,
    model_filter = () => true,
    models = Object.values(MODELS),
    point_style = {},
    show_model_labels = true,
    on_point_click,
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    x_prop: Label
    y_prop: Label
    models?: ModelData[]
    model_filter?: (model: ModelData) => boolean
    point_style?: PointStyle
    show_model_labels?: boolean
    on_point_click?: ComponentProps<typeof ScatterPlot>[`on_point_click`]
  } = $props()

  const label_path = (label: Label) => `${label?.path ?? ``}.${label?.key}`.replace(/^\./, ``)
  let x_path = $derived(label_path(x_prop))
  let y_path = $derived(label_path(y_prop))
  let x_label = $derived(x_prop?.label ?? x_prop?.key ?? `X`)
  let y_label = $derived(y_prop?.label ?? y_prop?.key ?? `Y`)

  let plot_data = $derived(
    models
      .filter(model_filter)
      .map((model) => ({
        x: get_nested_value(model, x_path),
        y: get_nested_value(model, y_path),
        metadata: { model_name: model.model_name, model_key: model.model_key },
        color: model?.color,
      }))
      .filter((pt) => pt.x != null && pt.y != null),
  )

  let series: DataSeries = $derived({
    x: plot_data.map((item) => Number(item.x)),
    y: plot_data.map((item) => Number(item.y)),
    markers: `points` as const,
    metadata: plot_data.map((item) => item.metadata),
    point_style: plot_data.map((item) => ({
      fill: item.color ?? `#4dabf7`,
      radius: 6,
      stroke: `white`,
      stroke_width: 0.5,
      ...point_style,
    })),
    // auto_placement repositions labels to avoid overlap, drawing leader lines back
    // to their points when displaced (see label_placement_config below)
    point_label: (show_model_labels ? plot_data : []).map((item) => ({
      text: item.metadata.model_name,
      font_size: `12px`,
      auto_placement: true,
    })),
  })
</script>

<ScatterPlot
  series={[series]}
  x_axis={{ label: x_label, format: x_prop?.format ?? `.1s`, range: [0, null] }}
  y_axis={{ label: y_label, format: y_prop?.format ?? `.3f`, range: [0, null] }}
  label_placement_config={{
    leader_line_threshold: 15,
    max_neighbors: { count: 3, radius: 40 },
  }}
  {on_point_click}
  {...rest}
>
  {#snippet tooltip({ x_formatted, y_formatted, metadata })}
    <strong>{metadata?.model_name}</strong><br />
    {x_label}: {x_formatted}<br />
    {y_label}: {y_formatted}<br />
  {/snippet}
</ScatterPlot>
