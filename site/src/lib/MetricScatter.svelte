<script lang="ts">
  import type { ModelData } from '$lib'
  import { calculate_days_ago, MODELS } from '$lib'
  import { get_nested_value } from '$lib/metrics'
  import type { Metric } from '$lib/types'
  import { ScatterPlot, type DataSeries, type PointStyle } from 'elementari'
  import { METADATA_COLS } from './labels'

  interface Props {
    x_prop: Metric
    y_prop: Metric
    models?: ModelData[]
    model_filter?: (model: ModelData) => boolean
    point_style?: PointStyle
    date_range?: [Date | null, Date | null]
    style?: string
    show_model_labels?: boolean | `auto-placement`
    [key: string]: unknown
  }
  let {
    x_prop,
    y_prop,
    model_filter = () => true,
    models = Object.values(MODELS),
    point_style = {},
    date_range = [null, null],
    style = ``,
    show_model_labels = true,
    ...rest
  }: Props = $props()

  // Add date range state for time series
  const now = new Date()
  const ms_per_day = 24 * 60 * 60 * 1000
  const n_days_ago = new Date(now.getTime() - 180 * ms_per_day)

  // Determine axis keys/paths
  let x_path = $derived(`${x_prop?.path ?? ``}.${x_prop?.key}`.replace(/^\./, ``))
  let y_path = $derived(`${y_prop?.path ?? ``}.${y_prop?.key}`.replace(/^\./, ``))
  // Determine labels
  let x_label = $derived(x_prop?.svg_label ?? x_prop?.label ?? x_prop?.key ?? `X`)
  let y_label = $derived(y_prop?.svg_label ?? y_prop?.label ?? y_prop?.key ?? `Y`)

  // Apply date range filter if needed
  function date_filter(model: ModelData): boolean {
    if (
      x_path !== METADATA_COLS.date_added.key ||
      (date_range[0] === null && date_range[1] === null)
    )
      return true

    const model_date = new Date(model.date_added ?? 0)
    return (
      model_date >= (date_range[0] ?? n_days_ago) && model_date <= (date_range[1] ?? now)
    )
  }

  // prepare and filter data
  let plot_data = $derived(
    models
      .filter((model) => model_filter(model) && date_filter(model))
      .map((model) => {
        let x_val = get_nested_value(model, x_path)
        if (x_path.includes(`date`)) x_val = new Date(x_val as string)
        let y_val = get_nested_value(model, y_path)
        if (y_path.includes(`date`)) y_val = new Date(y_val as string)
        const metadata = {
          model_name: model.model_name,
          date_added: model.date_added,
          days_ago: calculate_days_ago(model.date_added),
        }
        return { x: x_val, y: y_val, metadata, color: model?.color }
      })
      .filter((pt) => pt.x !== undefined && pt.y !== undefined),
  )

  // Create plot series based on show_model_labels
  let series: DataSeries = $derived({
    x: plot_data.map((item) => item.x) as number[],
    y: plot_data.map((item) => item.y) as number[],
    metadata: plot_data.map((item) => item.metadata),
    point_style: plot_data.map((item) => ({
      fill: item.color ?? `#4dabf7`,
      radius: 6,
      stroke: `white`,
      stroke_width: 0.5,
      ...point_style,
    })),
    point_label: (show_model_labels ? plot_data : []).map((item) => ({
      text: item.metadata.model_name,
      offset_y: 2,
      offset_x: 10,
      font_size: `12px`,
      color: item.color ?? `black`,
      auto_placement: show_model_labels === `auto-placement`,
    })),
  })
</script>

<ScatterPlot
  series={[series]}
  markers="points"
  {x_label}
  {y_label}
  x_format={x_prop?.format ?? `.1s`}
  y_format={y_prop?.format ?? `.3f`}
  {style}
  {...rest}
>
  {#snippet tooltip({ x_formatted, y_formatted, metadata })}
    {x_label}: {x_formatted}
    {#if x_path === `date_added`}
      <small>({metadata?.days_ago} days ago)</small>{/if}<br />
    {y_label}: {y_formatted}<br />
  {/snippet}
</ScatterPlot>
