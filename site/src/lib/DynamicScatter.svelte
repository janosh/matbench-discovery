<script lang="ts">
  import type { DiscoverySet, ModelData } from '$lib'
  import { calculate_days_ago } from '$lib'
  import { DEFAULT_CPS_CONFIG } from '$lib/combined_perf_score'
  import { format_property_path, get_format } from '$lib/properties'
  import type { DataSeries } from 'elementari'
  import { ScatterPlot } from 'elementari'
  import Select from 'svelte-multiselect'

  interface Props {
    models: ModelData[]
    model_filter?: (model: ModelData) => boolean
    discovery_set?: DiscoverySet
    point_color?: string | null
    point_radius?: number
    show_model_labels?: boolean
    style?: string
    [key: string]: unknown
  }
  let {
    models,
    model_filter = () => true,
    discovery_set = `unique_prototypes`,
    point_color = null,
    point_radius = 9,
    show_model_labels = true,
    style = ``,
    ...rest
  }: Props = $props()

  // Define property paths for available axes options
  const property_paths = [
    `CPS`,
    // Phonon metrics
    `phonons.kappa_103.κ_SRME`,
    // Geo opt metrics
    `geo_opt.symprec=1e-2.rmsd`,
    // Discovery metrics
    ...[`F1`, `Precision`, `Recall`, `Accuracy`, `DAF`, `R2`, `MAE`, `RMSE`].map(
      (key) => `discovery.${discovery_set}.${key}`,
    ),
    `model_params`,
    `date_added`,
    `n_training_materials`,
    `n_training_structures`,
    // Hyperparameter metrics for GNN and UIP models
    `hyperparams.graph_construction_radius`,
    `hyperparams.max_force`,
    `hyperparams.max_steps`,
    `hyperparams.batch_size`,
    `hyperparams.epochs`,
    `hyperparams.n_layers`,
  ]

  // Create lookup for property paths to their formatted labels
  const property_labels: Record<string, string> = {}
  property_paths.forEach((path) => {
    property_labels[path] = format_property_path(path)
  })

  // State for axis selection and log scale toggles
  let axes = $state({
    x: [`date_added`],
    y: [`CPS`],
    color: [`model_params`],
  })
  let log = $state({ x: false, y: false, color: true })

  // Format property paths for display with ">" separators
  let x_label = $derived(format_property_path(axes.x?.[0]))
  let y_label = $derived(format_property_path(axes.y?.[0]))
  let color_label = $derived(format_property_path(axes.color?.[0]))

  type MetricRecord = Record<
    string,
    Record<string, number | string | boolean | null | undefined>
  >

  // Extract property value from model with proper typing
  function get_nested_property(
    model: ModelData,
    property_path: string,
  ): number | undefined {
    // Handle special cases
    if (property_path === `model_params`)
      return typeof model.model_params === `number` ? model.model_params : undefined
    if (property_path === `n_estimators`)
      return typeof model.n_estimators === `number` ? model.n_estimators : undefined
    if (property_path === `date_added`)
      return model.date_added ? new Date(model.date_added).getTime() : undefined
    if (property_path === `n_training_materials`) return model.n_training_materials
    if (property_path === `n_training_structures`) return model.n_training_structures

    try {
      const parts = property_path.split(`.`)

      // Handle metrics pattern: category.set.metric_name
      if (parts.length === 3) {
        const [category, set_name, metric_name] = parts
        const metrics_obj = model.metrics?.[
          category as keyof typeof model.metrics
        ] as MetricRecord
        return typeof metrics_obj?.[set_name]?.[metric_name] === `number`
          ? (metrics_obj[set_name][metric_name] as number)
          : undefined
      }

      // Handle hyperparams pattern: hyperparams.param_name
      if (parts.length === 2 && parts[0] === `hyperparams`) {
        return typeof model.hyperparams?.[parts[1]] === `number`
          ? (model.hyperparams[parts[1]] as number)
          : undefined
      }

      // Generic fallback
      let value: unknown = model
      for (const part of parts) {
        if (value == null) return undefined
        value = value[part as keyof typeof value]
      }
      return typeof value === `number` ? value : undefined
    } catch {
      return undefined
    }
  }

  // Calculate counts for each property path across all models
  let model_counts_by_prop = $derived(
    property_paths.reduce(
      (acc, path) => {
        acc[path] = models.filter(
          (model) => get_nested_property(model, path) !== undefined,
        ).length
        return acc
      },
      {} as Record<string, number>,
    ),
  )

  // Prepare and filter data for display
  let filtered_models = $derived(
    models.filter((model) => {
      const x_val = get_nested_property(model, axes.x?.[0])
      const y_val = get_nested_property(model, axes.y?.[0])
      return (
        x_val !== undefined && // Filter models that have values for both x and y
        y_val !== undefined &&
        // AND color_value if point_color is null (meaning color scale is active)
        (point_color !== null ||
          get_nested_property(model, axes.color?.[0]) !== undefined) &&
        model_filter(model)
      )
    }),
  )

  let models_to_show = $derived(
    filtered_models.map((model) => {
      // For each property (x, y, color), get the value based on the selected property path
      const x = get_nested_property(model, axes.x?.[0])
      const y = get_nested_property(model, axes.y?.[0])
      const color_value = get_nested_property(model, axes.color?.[0])

      // Prepare metadata for display in tooltip
      const metadata = {
        model_name: model.model_name,
        date_added: model.date_added,
        days_ago: calculate_days_ago(model.date_added),
      }
      return { model, x: x ?? 0, y: y ?? 0, color_value, metadata }
    }),
  )

  // Check if any values are negative (for log scale controls)
  let x_has_neg_vals = $derived(models_to_show.some((item) => item.x <= 0))
  let y_has_neg_vals = $derived(models_to_show.some((item) => item.y <= 0))
  let color_has_neg_vals = $derived(
    models_to_show.some((item) => (item?.color_value ?? 1) <= 0),
  )
  // Auto-revert to linear scale when data contains negative or zero values
  $effect(() => {
    if (x_has_neg_vals) log.x = false
    if (y_has_neg_vals) log.y = false
    if (color_has_neg_vals) log.color = false
  })

  // Compute formats based on data
  let x_format = $derived(get_format(models_to_show.map((item) => item.x)))
  let y_format = $derived(get_format(models_to_show.map((item) => item.y)))

  // Set axis domains based on DEFAULT_CPS_CONFIG ranges if applicable
  function get_metric_range(prop_path: string | undefined): [number, number] | undefined {
    if (!prop_path) return
    const { label, range, parts } = DEFAULT_CPS_CONFIG

    if (prop_path === label && range) return range // selected prop is CPS

    // Otherwise, check if it's one of the parts
    const part = Object.values(parts).find((p) => p.path === prop_path)
    return part?.range
  }

  // Update series when dependencies change
  let series: DataSeries[] = $derived.by(() => {
    // Base styling for points
    const point_styles = models_to_show.map((item) => ({
      fill: point_color ?? item.model.color ?? `#4dabf7`,
      radius: point_radius,
      stroke: `white`,
      stroke_width: 0.5,
    }))

    let base_series: DataSeries = {
      x: models_to_show.map((item) => item.x),
      y: models_to_show.map((item) => item.y),
      point_style: point_styles,
      metadata: models_to_show.map((item) => item.metadata),
      color_values:
        point_color === null
          ? models_to_show
              .map((item) => item.color_value)
              .filter((v): v is number => v !== undefined)
          : undefined,
    }

    if (!show_model_labels) return [base_series]

    base_series.point_label = models_to_show.map((item) => ({
      text: item.metadata.model_name,
      offset_y: -5,
      offset_x: 5,
      font_size: `14px`,
      color: `black`,
      auto_placement: true,
    }))

    return [base_series]
  })
</script>

<div class="controls-grid">
  <!-- prettier-ignore -->
  {#each [
    { id: `x`, label: `X Axis`, selected_prop: axes.x?.[0], log_state: log.x, prop_value: axes.x?.[0], log_disabled: x_has_neg_vals },
    { id: `y`, label: `Y Axis`, selected_prop: axes.y?.[0], log_state: log.y, prop_value: axes.y?.[0], log_disabled: y_has_neg_vals },
    { id: `color`, label: `Color`, selected_prop: axes.color?.[0], log_state: log.color, prop_value: axes.color?.[0], log_disabled: color_has_neg_vals },
  ] as control (control.id)}
    <label for={control.id}>{control.label}</label>
    <Select
      id={control.id}
      bind:selected={axes[control.id as keyof typeof axes]}
      placeholder={`Select ${control.label}`}
      options={property_paths}
      maxSelect={1}
      minSelect={1}
      style="width: 100%; max-width: none; margin: 0;"
      liSelectedStyle="font-size: 16px;"
      ulSelectedStyle="padding: 0;"
      let:option
    >
      {@html property_labels[option]}
      <span style="font-size: smaller; color: gray; margin-left: 0.5em;">
        ({model_counts_by_prop[option]} models)
      </span>
    </Select>
    <label class:invisible={control.log_disabled || control.prop_value === `date_added`}>
      <input
        type="checkbox"
        bind:checked={log[control.id as keyof typeof log]}
        disabled={control.log_disabled}
      />
      Log {control.label.split(` `)[0].toLowerCase()} scale
    </label>
  {/each}
</div>

<div class="full-bleed-1400" style="height: 600px; margin-block: 1em;">
  {#key [log.color, axes.color?.[0]]}
    <ScatterPlot
      {series}
      {x_label}
      {y_label}
      {x_format}
      x_lim={get_metric_range(axes.x?.[0])}
      y_lim={get_metric_range(axes.y?.[0])}
      {y_format}
      markers="points"
      {style}
      x_scale_type={log.x ? `log` : `linear`}
      y_scale_type={log.y ? `log` : `linear`}
      x_label_shift={{ x: 0, y: -60 }}
      y_label_shift={{ x: 0, y: -20 }}
      color_scale_type={log.color ? `log` : `linear`}
      color_scheme="viridis"
      color_bar={{ label: color_label, label_side: `top` }}
      label_placement_config={{ link_strength: 10.0, link_distance: 5 }}
      {...rest}
    >
      {#snippet tooltip({ x_formatted, y_formatted, metadata })}
        <div
          style="white-space: nowrap; font-size: 14px; margin-top: 10px; line-height: 1;"
        >
          <strong>{metadata?.model_name}</strong><br />
          {@html x_label}: {x_formatted}
          {#if axes.x?.[0] === `date_added` && metadata?.days_ago}
            <small>({metadata.days_ago} days ago)</small>{/if}<br />
          {@html y_label.split(` > `).pop()}: {y_formatted}<br />
          {#if ![`model_params`, `date_added`].includes(axes.color?.[0]) && metadata?.model_name}
            {@html color_label.split(` > `).pop()}:
            {models_to_show
              .find((m) => m.metadata.model_name === metadata.model_name)
              ?.color_value?.toFixed(2) ?? `N/A`}<br />
          {/if}
        </div>
      {/snippet}
    </ScatterPlot>
  {/key}
</div>

<style>
  .invisible {
    visibility: hidden;
  }
  .controls-grid {
    display: grid;
    /* Add column for labels: Label Select Checkbox */
    grid-template-columns: auto 1fr auto;
    gap: 1em;
    align-items: center;
  }
  .controls-grid label {
    color: gray;
  }
  .controls-grid label[for] {
    text-align: right;
  }
  input[type='checkbox'] {
    transform: scale(1.2);
  }
</style>
