<script lang="ts">
  import type { CombinedMetricConfig, DiscoverySet, ModelData } from '$lib'
  import { calculate_combined_score } from '$lib/metrics'
  import { ScatterPlot, type DataSeries, type Point } from 'elementari'

  export type ModelProperty =
    | `model_params`
    | `training_cost`
    | `inference_cost`
    | `n_estimators`
    | `date_added`

  export type ModelMetric = `F1` | `RMSD` | `kappa_SRME` | `CPS` | string

  interface Props {
    models: ModelData[]
    config?: CombinedMetricConfig
    hovered?: boolean
    model_filter?: (model: ModelData) => boolean
    x_lim?: [number | null, number | null]
    y_lim?: [number | null, number | null] | null
    x_property?: ModelProperty
    y_metric?: ModelMetric
    metric?: string // Direct metric path for dotted notation (for backward compatibility)
    discovery_set?: DiscoverySet
    point_color?: string | null
    point_radius?: number
    x_label?: string
    y_label?: string
    tooltip_point?: Point | null
    date_range?: [Date | null, Date | null]
    style?: string
    x_format?: string
    [key: string]: unknown
  }

  let {
    models,
    config = undefined,
    hovered = $bindable(false),
    model_filter = () => true,
    y_lim = [0, 1],
    x_property = `model_params`,
    y_metric = `CPS`,
    metric = ``, // For dotted path access (MetricScatter compatibility)
    discovery_set = `unique_prototypes`,
    point_color = null,
    point_radius = 5,
    x_label,
    y_label,
    tooltip_point = $bindable(null),
    date_range = [null, null],
    style = ``,
    x_format = `.1s`,
    ...rest
  }: Props = $props()

  // Add date range state for time series
  const now = new Date()
  const ms_per_day = 24 * 60 * 60 * 1000
  const n_days_ago = new Date(now.getTime() - 180 * ms_per_day)

  // Define label lookups
  const property_labels: Record<string, string> = {
    model_params: `Model Parameters`,
    training_cost: `Training Cost (GPU hours)`,
    inference_cost: `Inference Cost`,
    n_estimators: `Number of Estimators`,
    date_added: `Date Added`,
  }

  const metric_labels: Record<string, string> = {
    F1: `F1 Score`,
    RMSD: `RMSD`,
    kappa_SRME: `κ SRME`,
    CPS: `Combined Performance Score`,
    // Default labels for dotted paths
    'phonons.kappa_103.κ_SRME': `κ SRME`,
    'discovery.unique_prototypes.F1': `F1 Score`,
    combined_performance_score: `Combined Performance Score`,
  }

  // If metric is provided directly (MetricScatter style), use it for y_metric
  let effective_y_metric = $derived(metric || y_metric)

  // Set default colors based on metric
  const default_color_map: Record<string, string> = {
    'phonons.kappa_103.κ_SRME': `#ff6b6b`,
    'discovery.unique_prototypes.F1': `#ffbb54`,
    combined_performance_score: `#a78bfa`,
    F1: `#ffbb54`,
    RMSD: `#4dabf7`,
    kappa_SRME: `#ff6b6b`,
    CPS: `#a78bfa`,
  }

  // Default colors and formats
  let default_color = $derived(default_color_map[effective_y_metric] ?? `#4dabf7`)
  let effective_color = $derived(point_color !== null ? point_color : default_color)
  let effective_x_format = $derived(x_property === `date_added` ? `%b %y` : x_format)

  // Determine labels
  let actual_x_label = $derived(x_label ?? property_labels[x_property] ?? `X`)
  let actual_y_label = $derived(
    y_label ?? metric_labels[effective_y_metric] ?? effective_y_metric ?? `Y`,
  )

  // Extract property values
  function get_property_value(
    model: ModelData,
    property: ModelProperty,
  ): number | undefined {
    if (property === `model_params`)
      return typeof model.model_params === `number` ? model.model_params : undefined
    if (property === `n_estimators`)
      return typeof model.n_estimators === `number` ? model.n_estimators : undefined
    if (property === `date_added`)
      return model.date_added ? new Date(model.date_added).getTime() : undefined

    if (
      property === `training_cost` &&
      typeof model.training_cost === `object` &&
      model.training_cost !== null
    ) {
      let total_gpu_hours = 0
      for (const [key, value] of Object.entries(model.training_cost)) {
        if (key.includes(`GPU`) && typeof value === `object` && value.hours) {
          total_gpu_hours += value.amount * value.hours
        }
      }
      return total_gpu_hours > 0 ? total_gpu_hours : undefined
    }

    return undefined
  }

  // Extract metric values, handling both enum values and dotted path strings
  function get_metric_value(
    model: ModelData,
    metric_key: ModelMetric,
  ): number | undefined {
    // If it's a standard enum metric, use the specific accessor logic
    if (metric_key === `F1`) return model.metrics?.discovery?.[discovery_set]?.F1

    if (
      metric_key === `RMSD` &&
      model.metrics?.geo_opt &&
      typeof model.metrics.geo_opt !== `string`
    ) {
      return model.metrics.geo_opt[`symprec=1e-5`]?.rmsd
    }

    if (
      metric_key === `kappa_SRME` &&
      model.metrics?.phonons &&
      typeof model.metrics.phonons !== `string`
    ) {
      return model.metrics.phonons.kappa_103?.κ_SRME !== undefined
        ? Number(model.metrics.phonons.kappa_103.κ_SRME)
        : undefined
    }

    if (metric_key === `CPS` && config) {
      const f1 = model.metrics?.discovery?.[discovery_set]?.F1
      const rmsd =
        model.metrics?.geo_opt && typeof model.metrics.geo_opt !== `string`
          ? model.metrics.geo_opt[`symprec=1e-5`]?.rmsd
          : undefined
      const kappa =
        model.metrics?.phonons && typeof model.metrics.phonons !== `string`
          ? model.metrics.phonons.kappa_103?.κ_SRME !== undefined
            ? Number(model.metrics.phonons.kappa_103.κ_SRME)
            : undefined
          : undefined

      const score = calculate_combined_score(f1, rmsd, kappa, config)
      return typeof score === `number` ? score : undefined
    }

    // If it's a dotted path string (from MetricScatter), use the nested accessor
    const nested_val = metric_key.split(`.`).reduce<unknown>((acc, key) => {
      if (acc && typeof acc === `object`) {
        return (acc as Record<string, unknown>)[key]
      }
      return undefined
    }, model.metrics)
    if (nested_val !== undefined) {
      return typeof nested_val === `number` ? nested_val : Number(nested_val)
    }

    return undefined
  }

  // Apply date range filter if needed
  function date_filter(model: ModelData): boolean {
    if (x_property !== `date_added` || (date_range[0] === null && date_range[1] === null))
      return true

    const model_date = new Date(model.date_added ?? 0)
    return (
      model_date >= (date_range[0] ?? n_days_ago) && model_date <= (date_range[1] ?? now)
    )
  }

  // Filter and prepare data
  let filtered_models = $derived(
    models
      .filter((model) => {
        const x_value = get_property_value(model, x_property)
        const y_value = get_metric_value(model, effective_y_metric)
        return (
          x_value !== undefined &&
          y_value !== undefined &&
          model_filter(model) &&
          date_filter(model)
        )
      })
      .sort((model1, model2) => {
        // For date-based plots, sort by date
        if (x_property === `date_added`) {
          const date_a = new Date(model1.date_added ?? 0)
          const date_b = new Date(model2.date_added ?? 0)
          return date_a.getTime() - date_b.getTime()
        }
        return 0 // No sorting for other properties
      }),
  )

  let plotable_models = $derived(
    filtered_models
      .map((model) => {
        const x_value = get_property_value(model, x_property)
        const y_value = get_metric_value(model, effective_y_metric)
        return {
          model,
          x: x_value !== undefined ? x_value : 0,
          y: y_value !== undefined ? y_value : 0,
          metadata: { model_name: model.model_name, date_added: model.date_added },
        }
      })
      .filter((item) => {
        // For properties that need positive values, filter out non-positive
        if (x_property !== `date_added`) {
          return item.x > 0 && item.y > 0
        }
        // For dates, just ensure y is valid
        return item.y > 0
      }),
  )

  // Find active model by matching tooltip point (for time series)
  let active_model = $derived(
    tooltip_point && x_property === `date_added`
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

  // Create point style
  let point_style = $derived({
    fill: effective_color,
    radius: point_radius,
    stroke: `white`,
    stroke_width: 0.5,
  })

  // Create plot series
  let series: DataSeries[] = $derived([
    {
      x: plotable_models.map((item) => item.x),
      y: plotable_models.map((item) => item.y),
      point_style,
      metadata: plotable_models.map((item) => item.metadata),
    },
  ])
</script>

<ScatterPlot
  {series}
  x_label={actual_x_label}
  y_label={actual_y_label}
  x_format={effective_x_format}
  y_format=".3f"
  bind:hovered
  bind:tooltip_point
  markers="points"
  {style}
  y_lim={y_lim ?? undefined}
  {...rest}
>
  {#snippet tooltip({ x_formatted, y_formatted, metadata })}
    {#if x_property === `date_added` && active_model}
      <div style="min-width: 10em;">
        <strong>{active_model.model_name}</strong>
        <br />
        {@html actual_y_label} = {@html y_formatted}
        <br />
        Added: {x_formatted}
      </div>
    {:else}
      <div style="white-space: nowrap;">
        <strong>{metadata?.model_name}</strong><br />
        {actual_x_label}: {x_formatted}<br />
        {actual_y_label}: {y_formatted}
      </div>
    {/if}
  {/snippet}
</ScatterPlot>
