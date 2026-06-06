<script lang="ts">
  import type { ModelData } from '$lib'
  import { MODELS } from '$lib'
  import { wide_legend } from '$lib/fig-helpers'
  import { get_nested_value, is_finite_num } from '$lib/metrics'
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

  // One series per model so matterviz's built-in legend can toggle models individually and
  // auto-assign a distinct marker shape per series (colored by each model's brand color)
  let series: DataSeries[] = $derived(
    models
      .filter(model_filter)
      .map((model) => ({
        model,
        x: get_nested_value(model, x_path),
        y: get_nested_value(model, y_path),
      }))
      .filter((pt) => is_finite_num(pt.x) && is_finite_num(pt.y))
      .map(({ model, x, y }) => ({
        x: [Number(x)],
        y: [Number(y)],
        label: model.model_name,
        markers: `points` as const,
        metadata: [{ model_name: model.model_name, model_key: model.model_key }],
        point_style: {
          fill: model.color ?? `#4dabf7`,
          radius: 6,
          stroke: `white`,
          stroke_width: 0.5,
          ...point_style,
        },
        // auto_placement repositions labels to avoid overlap, drawing leader lines back
        // to their points when displaced (see label_placement_config below)
        point_label: show_model_labels
          ? [{ text: model.model_name, font_size: `12px`, auto_placement: true }]
          : [],
      })),
  )
</script>

<ScatterPlot
  series={series}
  x_axis={{ label: x_label, format: x_prop?.format ?? `.1s`, range: [0, null] }}
  y_axis={{ label: y_label, format: y_prop?.format ?? `.3f`, range: [0, null] }}
  legend={wide_legend}
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
