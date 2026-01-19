<script lang="ts">
  import { goto } from '$app/navigation'
  import type { ModelData } from '$lib'
  import { calculate_days_ago } from '$lib'
  import { extent } from 'd3-array'
  import { DraggablePane, format_num, Icon } from 'matterviz'
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

  const date_key = METADATA_COLS.date_added.key
  const params_key = HYPERPARAMS.model_params.key

  const { model_params, graph_construction_radius, max_force, max_steps } =
    HYPERPARAMS
  const { batch_size, epochs, n_layers } = HYPERPARAMS
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

  // Color and Size are managed separately (not part of ScatterPlot's built-in axis selection)
  let color_key = $state(ALL_METRICS.F1.key)
  let size_prop = $state(HYPERPARAMS.model_params)

  // Derive color_prop from key
  let color_prop = $derived(options_by_key[color_key])

  // Derive full axes object for use in plot_data computation
  let axes = $derived({
    x: options_by_key[selected.x],
    y: options_by_key[selected.y],
    color_value: color_prop,
    size_value: size_prop,
  })

  let log = $state({ x: false, y: false, color: false, size: false })
  let is_fullscreen = $state(false)
  let show_extra_controls = $state(false)
  let container_el: HTMLDivElement | null = null

  let color_scheme: D3InterpolateName = $state(`interpolateViridis`)
  let x_ticks = $state(5)
  let y_ticks = $state(5)
  let x_grid = $state(true)
  let y_grid = $state(true)

  let size_multiplier = $state(1)
  let label_font_size = $state(14)
  let link_strength = $state(5)
  let link_distance = $state(5)
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

  // Convert options to format expected by ScatterPlot axis options
  let axis_options = $derived(
    options.map((prop) => ({
      key: prop.key,
      label: `${prop.label} (${model_counts_by_prop[prop.key]} models)`,
    })),
  )

  // Options for ColorBar property select
  let color_options = $derived(
    options.map((prop) => ({
      key: prop.key,
      label: `${prop.label} (${model_counts_by_prop[prop.key]} models)`,
      unit: prop.unit,
    })),
  )

  // Data loader for interactive axis selection (y2 not supported in this scatter plot)
  const data_loader = async (axis: `x` | `y` | `y2`, key: string) => {
    // Update the selected key
    if (axis === `x`) selected.x = key
    else if (axis === `y`) selected.y = key

    // Wait for reactive derivations to update
    await tick()

    // Return updated series
    return { series: [series], axis_label: options_by_key[key]?.label }
  }

  function is_num_or_date(val: unknown): boolean {
    if (typeof val === `number` && !isNaN(val)) return true
    if (val instanceof Date) return !isNaN(val.getTime())
    return false
  }

  let plot_data = $derived(
    models
      .filter(model_filter)
      .map((model) => {
        const x_val = get_label_value(model, axes.x)
        const y_val = get_label_value(model, axes.y)
        let color_value = get_label_value(model, axes.color_value)
        let size_value = get_label_value(model, axes.size_value)
        if (axes.size_value?.key === date_key) {
          const timestamp = new Date(String(size_value)).getTime()
          if (!isNaN(timestamp)) size_value = timestamp
        }

        const { model_name, date_added, color, model_key } = model
        const days_ago = calculate_days_ago(model.date_added)
        const metadata = { model_name, date_added, days_ago, model_key }
        return { x: x_val, y: y_val, color_value, size_value, metadata, color }
      })
      .filter((item) => {
        const x_valid = is_num_or_date(item.x)
        const y_valid = is_num_or_date(item.y)
        const color_valid = is_num_or_date(item.color_value)
        const size_valid = is_num_or_date(item.size_value)
        return x_valid && y_valid && color_valid && size_valid
      }),
  )

  let series = $derived({
    x: plot_data.map((item) => item.x as number),
    y: plot_data.map((item) => item.y as number),
    markers: `points` as const,
    point_style: plot_data.map((item) => ({ fill: point_color ?? item.color })),
    metadata: plot_data.map((item) => item.metadata),
    color_values: point_color === null
      ? (plot_data
        .map((item) => item.color_value)
        .filter((val): val is number | Date => val !== undefined)
        .map((val) => (val instanceof Date ? val.getTime() : val)) as number[])
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

<svelte:window
  onfullscreenchange={() => is_fullscreen = document.fullscreenElement === container_el}
/>

<div
  bind:this={container_el}
  class:fullscreen={is_fullscreen}
  class:bleed-1400={!is_fullscreen}
  style:height={is_fullscreen ? `100%` : `auto`}
  style="margin-block: 2em"
>
  <div class="controls-row">
    <label for="size-select">Marker Size</label>
    <Select
      {options}
      id="size-select"
      bind:value={size_prop}
      placeholder="Select Size"
      maxSelect={1}
      minSelect={1}
      style="flex: 1; max-width: 300px; margin: 0"
      liSelectedStyle="font-size: 14px;"
    >
      {#snippet children({ option: prop }: { option: typeof options[number] })}
        {@html format_property_path(
          `${prop.path ?? ``}.${prop.key}`.replace(/^\./, ``),
        )}
        <span style="font-size: smaller; color: gray; margin-left: 0.5em">
          ({model_counts_by_prop[prop.key]} models)
        </span>
      {/snippet}
    </Select>
  </div>

  <ScatterPlot
    style={`height: ${is_fullscreen ? `100%` : `600px`}`}
    series={[series]}
    x_axis={{
      label: axes.x?.label,
      format: axes.x?.format,
      range: axes.x?.range,
      scale_type: log.x ? `log` : `linear`,
      label_shift: { y: -50 },
      ticks: x_ticks,
      options: axis_options,
      selected_key: selected.x,
    }}
    y_axis={{
      label: axes.y?.label,
      format: axes.y?.format,
      range: axes.y?.range,
      scale_type: log.y ? `log` : `linear`,
      label_shift: {
        x: 50,
        y: [date_key, params_key].includes(axes.y?.key ?? ``) ? -40 : -10,
      },
      ticks: y_ticks,
      options: axis_options,
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
      property_options: color_options,
      selected_property_key: color_key,
      data_loader: async (key) => {
        color_key = key
        const prop = options_by_key[key]
        const values = models
          .filter(model_filter)
          .map((model) => get_label_value(model, prop))
          .filter((val): val is number => typeof val === `number` && !isNaN(val))
        const [min, max] = extent(values) as [number, number]
        const title = `${prop?.label}${
          prop?.better ? ` (${prop?.better}=better)` : ``
        }`
        return { range: [min ?? 0, max ?? 1], title }
      },
    }}
    label_placement_config={{
      link_strength,
      link_distance,
      link_distance_range: [min_link_distance, max_link_distance],
    }}
    point_events={{
      onclick: ({ point }) => goto(`/models/${point.metadata?.model_key ?? ``}`),
    }}
    {...rest}
    {data_loader}
  >
    <DraggablePane
      bind:show={show_extra_controls}
      toggle_props={{
        class: `settings-toggle`,
        style:
          `position: absolute; top: 15px; right: 4.3em; background-color: var(--btn-bg); border-radius: 50%; padding: 4pt;`,
      }}
      pane_props={{
        style: `border: 1px solid var(--border);`,
        'aria-hidden': !show_extra_controls,
      }}
    >
      <div style="display: grid; grid-template-columns: auto 1fr; gap: 8pt 1em">
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
            style="grid-column: 1/-1; display: flex; gap: 1em; flex-wrap: wrap"
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

        <label
          style="grid-column: 1/-1"
          title="Toggle visibility of model name labels on the scatter plot points"
        >
          <input type="checkbox" bind:checked={show_model_labels} /> Show Labels
        </label>
        <ColorScaleSelect
          bind:value={color_scheme}
          style="margin: 0; grid-column: 1/-1"
        />
        <label title="Toggle the visibility of vertical grid lines">
          <input type="checkbox" bind:checked={x_grid} /> X Grid
        </label>
        <label title="Set the approximate number of ticks on the X axis">Ticks:
          <input
            id="x-ticks"
            type="number"
            min="0"
            max="20"
            bind:value={x_ticks}
            style="width: 50px"
          /></label>

        <label title="Toggle the visibility of horizontal grid lines">
          <input type="checkbox" bind:checked={y_grid} /> Y Grid
        </label>
        <label title="Set the approximate number of ticks on the Y axis">Ticks:
          <input
            id="y-ticks"
            type="number"
            min="0"
            max="20"
            bind:value={y_ticks}
            style="width: 50px"
          />
        </label>
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
          style="grid-column: 2"
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
      </div>
    </DraggablePane>

    <button
      onclick={() => {
        if (document.fullscreenElement === container_el) {
          document.exitFullscreen()
          is_fullscreen = false
        } else {
          container_el?.requestFullscreen?.()
          is_fullscreen = true
        }
      }}
      aria-label="{is_fullscreen ? `Exit` : `Enter`} fullscreen"
      title="{is_fullscreen ? `Exit` : `Enter`} fullscreen"
    >
      <Icon icon="{is_fullscreen ? `Exit` : ``}Fullscreen" />
    </button>

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
  div.fullscreen {
    background: var(--page-bg);
    overflow: auto;
    padding: 2em;
    box-sizing: border-box;
    display: flex;
    flex-direction: column;
    gap: 1em;
  }
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
  button[title$='fullscreen'] {
    position: absolute;
    top: 1em;
    right: 3.3em;
    display: flex;
    padding: 8px;
    border-radius: 50%;
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
