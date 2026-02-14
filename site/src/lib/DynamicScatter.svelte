<script lang="ts">
  import { goto } from '$app/navigation'
  import type { ModelData } from '$lib'
  import { calculate_days_ago } from '$lib'
  import { extent } from 'd3-array'
  import { format_num } from 'matterviz'
  import type { D3InterpolateName } from 'matterviz/colors'
  import { ColorScaleSelect, ScatterPlot } from 'matterviz/plot'
  import type { ComponentProps } from 'svelte'
  import { tick } from 'svelte'
  import Select from 'svelte-multiselect'
  import {
    ALL_METRICS,
    format_property_path,
    HYPERPARAMS,
    METADATA_COLS,
  } from './labels'
  import { get_nested_value } from './metrics'
  import type { Label } from './types'

  // Build data access path from label, handling property vs key distinction
  function get_label_path(label: Label | undefined): string {
    const prop_name = label?.property ?? label?.key ?? ``
    return `${label?.path ?? ``}.${prop_name}`.replace(/^\./, ``)
  }

  // Get value from model using label's path, converting dates to timestamps
  function get_label_value(model: ModelData, label: Label | undefined): unknown {
    const path = get_label_path(label)
    let val = get_nested_value(model, path)
    if (path.includes(`date`)) val = new Date(val as string).getTime()
    return val
  }

  let {
    models,
    model_filter = () => true,
    point_color = null,
    show_model_labels = true,
    ...rest
  }: ComponentProps<typeof ScatterPlot> & {
    models: ModelData[]
    model_filter?: (model: ModelData) => boolean
    point_color?: string | null
    show_model_labels?: boolean
  } = $props()

  const {
    model_params,
    graph_construction_radius,
    max_force,
    max_steps,
    batch_size,
    epochs,
    n_layers,
  } = HYPERPARAMS
  const { date_added, n_training_materials, n_training_structures } = METADATA_COLS

  const options = [
    ...Object.values(ALL_METRICS),
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

  // Create lookup map from key to full property object
  const options_by_key = Object.fromEntries(options.map((opt) => [opt.key, opt]))

  // Track selected axis keys for x/y (used by ScatterPlot's built-in axis selection)
  let selected = $state({
    x: ALL_METRICS.κ_SRME.key,
    y: ALL_METRICS.CPS.key,
  })

  let color_key = $state(ALL_METRICS.F1.key)
  let size_prop = $state(HYPERPARAMS.model_params as typeof options[number])

  let axes = $derived({
    x: options_by_key[selected.x],
    y: options_by_key[selected.y],
    color_value: options_by_key[color_key],
    size_value: size_prop,
  })

  let log = $state({ x: false, y: false, color: false, size: false })

  let color_scheme: D3InterpolateName = $state(`interpolateViridis`)
  // Grid/tick settings are initial values — PlotControls handles the UI
  const [x_ticks, y_ticks, x_grid, y_grid] = [5, 5, true, true]

  let size_multiplier = $state(1)
  let label_font_size = $state(14)
  let link_strength = $state(5)
  let min_link_distance = $state(15)
  let max_link_distance = $state(20)

  let model_counts_by_prop = $derived(
    options.reduce(
      (acc, prop) => {
        const path = get_label_path(prop)
        acc[prop.key] = models.filter(
          (model) => get_nested_value(model, path) !== undefined,
        ).length
        return acc
      },
      {} as Record<string, number>,
    ),
  )

  // Options for axis selects and color bar property select
  let prop_options = $derived(
    options.map((prop) => ({
      key: prop.key,
      label: `${prop.label} (${model_counts_by_prop[prop.key]} models)`,
      unit: prop.unit,
    })),
  )

  // Data loader for interactive axis selection (y2 not supported in this scatter plot)
  const data_loader = async (axis: `x` | `y` | `y2`, key: string) => {
    if (axis === `x`) selected.x = key
    else if (axis === `y`) selected.y = key

    await tick()
    return { series: [series], axis_label: options_by_key[key]?.label }
  }

  // get_label_value already converts dates to timestamps, so only numbers reach here
  const is_finite_num = (val: unknown): val is number =>
    typeof val === `number` && isFinite(val)

  let plot_data = $derived(
    models
      .filter(model_filter)
      .map((model) => {
        const x_val = get_label_value(model, axes.x)
        const y_val = get_label_value(model, axes.y)
        const color_value = get_label_value(model, axes.color_value)
        const size_value = get_label_value(model, axes.size_value)

        const { model_name, date_added: model_date, color, model_key } = model
        const days_ago = calculate_days_ago(model.date_added)
        const metadata = { model_name, date_added: model_date, days_ago, model_key }
        return { x: x_val, y: y_val, color_value, size_value, metadata, color }
      })
      .filter((item) => {
        const required = [item.x, item.y, item.size_value]
        // only require finite color_value when no fixed point_color is set
        if (point_color === null) required.push(item.color_value)
        return required.every(is_finite_num)
      }),
  )

  let series = $derived({
    x: plot_data.map((item) => item.x as number),
    y: plot_data.map((item) => item.y as number),
    markers: `points` as const,
    point_style: plot_data.map((item) => ({ fill: point_color ?? item.color })),
    metadata: plot_data.map((item) => item.metadata),
    color_values: point_color === null
      ? plot_data.map((item) => item.color_value as number)
      : undefined,
    size_values: axes.size_value
      ? plot_data.map((item) => item.size_value as number)
      : undefined,
    point_label: (show_model_labels ? plot_data : []).map((item) => ({
      text: item.metadata.model_name,
      offset_y: 0,
      offset_x: 10,
      font_size: `${label_font_size}px`,
      color: `black`,
      auto_placement: true,
    })),
  })
</script>

<div class="bleed-1400" style="margin-block: 2em">
  <div class="controls-row">
    <label for="size-select">Marker Size</label>
    <Select
      {options}
      id="size-select"
      bind:value={size_prop}
      maxSelect={1}
      minSelect={1}
      style="flex: 1; max-width: 300px; margin: 0"
      liSelectedStyle="font-size: 14px;"
    >
      {#snippet children({ option: prop }: { option: typeof options[number] })}
        {@html format_property_path(get_label_path(prop))}
        <span style="font-size: smaller; color: gray; margin-left: 0.5em">
          ({model_counts_by_prop[prop.key]} models)
        </span>
      {/snippet}
    </Select>
  </div>

  <ScatterPlot
    style="height: 600px"
    series={[series]}
    x_axis={{
      label: axes.x?.label,
      format: axes.x?.format,
      range: axes.x?.range,
      scale_type: log.x ? `log` : `linear`,
      label_shift: { y: -50 },
      ticks: x_ticks,
      options: prop_options,
      selected_key: selected.x,
    }}
    y_axis={{
      label: axes.y?.label,
      format: axes.y?.format,
      range: axes.y?.range,
      scale_type: log.y ? `log` : `linear`,
      label_shift: {
        x: 50,
        y: [`date_added`, `model_params`].includes(axes.y?.key ?? ``) ? -40 : -10,
      },
      ticks: y_ticks,
      options: prop_options,
      selected_key: selected.y,
    }}
    display={{ x_grid, y_grid }}
    color_scale={{ scheme: color_scheme, type: log.color ? `log` : `linear` }}
    size_scale={{
      radius_range: [5 * size_multiplier, 20 * size_multiplier],
      type: log.size ? `log` : `linear`,
    }}
    color_bar={{
      title: `${axes.color_value?.label}${
        axes.color_value?.better ? ` (${axes.color_value?.better}=better)` : ``
      }`,
      margin: { t: 30, l: 80, b: 80, r: 50 },
      tick_format: axes.color_value?.format,
      property_options: prop_options,
      selected_property_key: color_key,
      data_loader: async (key) => {
        color_key = key
        const prop = options_by_key[key]
        const values = models
          .filter(model_filter)
          .map((model) => get_label_value(model, prop))
          .filter((val): val is number => typeof val === `number` && isFinite(val))
        const [min, max] = extent(values) as [number, number]
        const title = `${prop?.label}${
          prop?.better ? ` (${prop?.better}=better)` : ``
        }`
        return { range: [min ?? 0, max ?? 1], title }
      },
    }}
    label_placement_config={{
      link_strength,
      link_distance_range: [min_link_distance, max_link_distance],
    }}
    point_events={{
      onclick: ({ point }) => goto(`/models/${point.metadata?.model_key ?? ``}`),
    }}
    {...rest}
    {data_loader}
  >
    {#snippet controls_extra()}
      <!-- Log scale toggles - {#if true} creates scope for {@const} declarations -->
      {#if true}
        {@const x_extent = extent(plot_data, (d) => d.x as number)}
        {@const y_extent = extent(plot_data, (d) => d.y as number)}
        {@const color_extent = extent(plot_data, (d) => d.color_value as number)}
        {@const size_extent = extent(plot_data, (d) => d.size_value as number)}
        {@const can_log_x = x_extent[0] !== undefined && x_extent[0] > 0 &&
        100 * x_extent[0] <= (x_extent[1] ?? 0)}
        {@const can_log_y = y_extent[0] !== undefined && y_extent[0] > 0 &&
        100 * y_extent[0] <= (y_extent[1] ?? 0)}
        {@const can_log_color = color_extent[0] !== undefined &&
        color_extent[0] > 0 &&
        100 * color_extent[0] <= (color_extent[1] ?? 0)}
        {@const can_log_size = size_extent[0] !== undefined && size_extent[0] > 0 &&
        100 * size_extent[0] <= (size_extent[1] ?? 0)}
        <div
          class="log-toggles"
          style="display: flex; gap: 1em; flex-wrap: wrap"
        >
          <label style:visibility={can_log_x ? `visible` : `hidden`}>
            <input type="checkbox" bind:checked={log.x} disabled={!can_log_x} /> Log X
          </label>
          <label style:visibility={can_log_y ? `visible` : `hidden`}>
            <input type="checkbox" bind:checked={log.y} disabled={!can_log_y} /> Log Y
          </label>
          <label style:visibility={can_log_color ? `visible` : `hidden`}>
            <input type="checkbox" bind:checked={log.color} disabled={!can_log_color} />
            Log Color
          </label>
          <label style:visibility={can_log_size ? `visible` : `hidden`}>
            <input type="checkbox" bind:checked={log.size} disabled={!can_log_size} />
            Log Size
          </label>
        </div>
      {/if}

      <label title="Toggle visibility of model name labels on the scatter plot points">
        <input type="checkbox" bind:checked={show_model_labels} /> Show Labels
      </label>
      <ColorScaleSelect
        bind:value={color_scheme}
        style="margin: 0"
      />
      <label
        for="size-multiplier"
        title="Adjust the base size of all points on the scatter plot (multiplier for radius)"
      >Point Size</label>
      <input
        id="size-multiplier"
        type="range"
        min="0.1"
        max="5"
        step="0.1"
        bind:value={size_multiplier}
      />
      <label
        for="label-font-size"
        title="Adjust the font size of the model name labels (in pixels)"
      >Label Size</label>
      <input
        id="label-font-size"
        type="range"
        min="8"
        max="24"
        step="1"
        bind:value={label_font_size}
      />
      <label
        title="Configure the distance range and strength of the links connecting labels to their points"
        for="min-link-distance"
      >Label Link</label>
      <div class="combined-link-controls">
        <input
          id="min-link-distance"
          type="number"
          min="0"
          max="100"
          bind:value={min_link_distance}
          title="Minimum distance"
        />
        <span>-</span>
        <input
          id="max-link-distance"
          type="number"
          min="0"
          max="100"
          bind:value={max_link_distance}
          title="Maximum distance"
        />
        <input
          id="link-strength"
          type="range"
          min="0.1"
          max="10"
          step="0.1"
          bind:value={link_strength}
          title="Strength (higher = stronger pull)"
          style="flex: 1"
        />
      </div>
    {/snippet}

    {#snippet tooltip({ x_formatted, y_formatted, metadata })}
      {#if metadata}
        {@const point = plot_data.find((m) => m.metadata.model_name === metadata.model_name)}
        <strong>{metadata.model_name}</strong><br />
        {@html axes.x?.label}: {x_formatted}
        {#if axes.x?.key === `date_added` && metadata.days_ago}
          <small>({metadata.days_ago} days ago)</small>{/if}<br />
        {@html axes.y?.label}: {y_formatted}<br />
        {#if ![`model_params`, `date_added`].includes(axes.color_value?.key ?? ``) &&
        point?.color_value !== undefined}
          {@html axes.color_value?.label}:
          {format_num(point.color_value as number)}<br />
        {/if}
        {#if axes.size_value && point?.size_value !== undefined}
          {@html axes.size_value.label}:
          {format_num(point.size_value as number)}<br />
        {/if}
      {/if}
    {/snippet}
  </ScatterPlot>
</div>

<style>
  div.controls-row {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 1ex 2em;
    margin: 0 3em 1em;
    justify-content: center;
  }
  div.controls-row label {
    font-weight: 500;
  }
  div.combined-link-controls {
    display: flex;
    gap: 0.5em;
    align-items: center;
  }
  div.combined-link-controls input[type='number'] {
    width: 50px;
  }
</style>
