<script lang="ts">
  import type { Metric, ModelData } from '$lib'
  import { calculate_days_ago } from '$lib'
  import type { DataSeries } from 'elementari'
  import { ScatterPlot } from 'elementari'
  import Select from 'svelte-multiselect'
  import { format_property_path, HYPERPARAMS, METADATA_COLS, METRICS } from './labels'
  import { get_nested_value } from './metrics'

  interface Props {
    models: ModelData[]
    model_filter?: (model: ModelData) => boolean
    point_color?: string | null
    point_radius?: number
    show_model_labels?: boolean
    style?: string
    [key: string]: unknown
  }
  let {
    models,
    model_filter = () => true,
    point_color = null,
    point_radius = 9,
    show_model_labels = true,
    style = ``,
    ...rest
  }: Props = $props()

  const date_key = METADATA_COLS.date_added.key
  const params_key = METADATA_COLS.model_params.key

  // State for axis selection and log scale toggles
  let axes = $state({
    x: METRICS.κ_SRME,
    y: METRICS.CPS,
    color: METRICS.F1,
  })
  let log = $state({ x: false, y: false, color: true })

  const { model_params, date_added, n_training_materials, n_training_structures } =
    METADATA_COLS
  const {
    graph_construction_radius,
    max_force,
    max_steps,
    batch_size,
    epochs,
    n_layers,
  } = HYPERPARAMS
  const options = [
    ...Object.values(METRICS),
    model_params,
    date_added,
    n_training_materials,
    n_training_structures,
    graph_construction_radius,
    max_force,
    max_steps,
    batch_size,
    epochs,
    n_layers,
  ]
  // Calculate counts for each property path across all models
  let model_counts_by_prop = $derived(
    options.reduce(
      (acc, prop) => {
        const path = `${prop.path ?? ``}.${prop.key}`.replace(/^\./, ``)
        acc[prop.key] = models.filter(
          (model) => get_nested_value(model, path) !== undefined,
        ).length
        return acc
      },
      {} as Record<string, number>,
    ),
  )

  function is_num_or_date(val: unknown): boolean {
    return (
      (typeof val === `number` && !isNaN(val)) || val instanceof Date // warning: this is true for invalid dates
    )
  }

  let plot_data = $derived(
    models
      .filter(model_filter)
      .map((model) => {
        // For each property (x, y, color), get the value based on the selected property path
        let x_path = `${axes.x.path ?? ``}.${axes.x.key}`
        if (x_path.startsWith(`.`)) x_path = x_path.slice(1)
        let x_val = get_nested_value(model, x_path)
        if (x_path.includes(`date`)) x_val = new Date(x_val as string)

        let y_path = `${axes.y.path ?? ``}.${axes.y.key}`
        if (y_path.startsWith(`.`)) y_path = y_path.slice(1)
        let y_val = get_nested_value(model, y_path)
        if (y_path.includes(`date`)) y_val = new Date(y_val as string)
        let color_path = `${axes.color.path ?? ``}.${axes.color.key}`
        if (color_path.startsWith(`.`)) color_path = color_path.slice(1)
        let color_value = get_nested_value(model, color_path)
        if (color_path.includes(`date`)) color_value = new Date(color_value as string)
        // Prepare metadata for display in tooltip
        const { model_name, date_added, color } = model
        const metadata = {
          model_name,
          date_added,
          days_ago: calculate_days_ago(model.date_added),
        }
        return { x: x_val, y: y_val, color_value, metadata, color }
      })
      .filter((item) => {
        const x_valid = is_num_or_date(item.x)
        const y_valid = is_num_or_date(item.y)
        const color_valid = is_num_or_date(item.color_value)
        return x_valid && y_valid && color_valid
      }),
  )

  // Check if any values are negative (for log scale controls)
  let x_has_neg_vals = $derived(plot_data.some((item) => item.x <= 0))
  let y_has_neg_vals = $derived(plot_data.some((item) => item.y <= 0))
  let color_has_neg_vals = $derived(
    plot_data.some((item) => (item?.color_value ?? 1) <= 0),
  )
  // Auto-revert to linear scale when data contains negative or zero values
  $effect(() => {
    if (x_has_neg_vals) log.x = false
    if (y_has_neg_vals) log.y = false
    if (color_has_neg_vals) log.color = false
  })

  // Update series when dependencies change
  let series: DataSeries[] = $derived.by(() => {
    // Base styling for points
    const point_styles = plot_data.map((item) => ({
      fill: point_color ?? item.color ?? `#4dabf7`,
      radius: point_radius,
      stroke: `white`,
      stroke_width: 0.5,
    }))

    let base_series: DataSeries = {
      x: plot_data.map((item) => item.x),
      y: plot_data.map((item) => item.y),
      point_style: point_styles,
      metadata: plot_data.map((item) => item.metadata),
      color_values:
        point_color === null
          ? plot_data
              .map((item) => item.color_value)
              .filter((v): v is number => v !== undefined)
          : undefined,
    }

    if (!show_model_labels) return [base_series]

    base_series.point_label = plot_data.map((item) => ({
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
    { id: `x`, label: `X Axis`, log_state: log.x, log_disabled: x_has_neg_vals },
    { id: `y`, label: `Y Axis`, log_state: log.y, log_disabled: y_has_neg_vals },
    { id: `color`, label: `Color`, log_state: log.color, log_disabled: color_has_neg_vals },
  ] as control (control.id)}
    <label for={control.id}>{control.label}</label>
    <Select
      id={control.id}
      selected={[axes[control.id as keyof typeof axes]]}
      bind:value={axes[control.id as keyof typeof axes]}
      placeholder="Select {control.label}"
      {options}
      maxSelect={1}
      minSelect={1}
      style="width: 100%; max-width: none; margin: 0;"
      liSelectedStyle="font-size: 16px;"
      ulSelectedStyle="padding: 0;"
      let:option
      --sms-selected-bg="none"
      --sms-border="1px solid rgba(255, 255, 255, 0.15)"
    >
      {@const prop = option as unknown as Metric}
      {@html format_property_path(`${prop.path ?? ``}.${prop.key}`.replace(/^\./, ``))}
      <span style="font-size: smaller; color: gray; margin-left: 0.5em;">
        ({model_counts_by_prop[prop.key]} models)
      </span>
    </Select>
    <label class:invisible={control.log_disabled || axes[control.id as keyof typeof axes]?.key === `date_added`}>
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
  {#key [log.color, axes.color.key]}
    <!-- TODO fix x_lim and y_lim to use metric ranges-->
    <ScatterPlot
      {series}
      x_label={axes.x.label}
      y_label={axes.y.label}
      x_lim={axes.x.range}
      y_lim={axes.y.range}
      x_format={axes.x.format}
      y_format={axes.y.format}
      markers="points"
      {style}
      x_scale_type={log.x ? `log` : `linear`}
      y_scale_type={log.y ? `log` : `linear`}
      x_label_shift={{ y: -60 }}
      y_label_shift={{ y: [date_key, params_key].includes(axes.y.key ?? ``) ? -40 : -10 }}
      color_scale_type={log.color ? `log` : `linear`}
      color_scheme="viridis"
      color_bar={{ label: axes.color.label, label_side: `top`, margin: 30 }}
      label_placement_config={{ link_strength: 2, link_distance: 1 }}
      {...rest}
    >
      {#snippet tooltip({ x_formatted, y_formatted, metadata })}
        <strong>{metadata?.model_name}</strong><br />
        {@html axes.x.label}: {x_formatted}
        {#if axes.x.key === `date_added` && metadata?.days_ago}
          <small>({metadata.days_ago} days ago)</small>{/if}<br />
        {@html axes.y.label}: {y_formatted}<br />
        {#if ![`model_params`, `date_added`].includes(axes.color.key ?? ``) && metadata?.model_name}
          {@html axes.color.label}:
          {plot_data
            .find((m) => m.metadata.model_name === metadata.model_name)
            ?.color_value?.toFixed(2) ?? `N/A`}<br />
        {/if}
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
