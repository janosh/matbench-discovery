<script lang="ts">
  import type { Metric, ModelData } from '$lib'
  import { calculate_days_ago } from '$lib'
  import { extent } from 'd3-array'
  import type { DataSeries } from 'elementari'
  import { ColorScaleSelect, pretty_num, ScatterPlot } from 'elementari'
  import type { D3ColorSchemeName } from 'elementari/colors'
  import Select from 'svelte-multiselect'
  import { click_outside } from 'svelte-zoo'
  import {
    ALL_METRICS,
    format_property_path,
    HYPERPARAMS,
    METADATA_COLS,
    to_title,
  } from './labels'
  import { get_nested_value } from './metrics'

  interface Props {
    models: ModelData[]
    model_filter?: (model: ModelData) => boolean
    point_color?: string | null
    show_model_labels?: boolean
    style?: string
    [key: string]: unknown
  }
  let {
    models,
    model_filter = () => true,
    point_color = null,
    show_model_labels = true,
    style = ``,
    ...rest
  }: Props = $props()

  const date_key = METADATA_COLS.date_added.key
  const params_key = HYPERPARAMS.model_params.key

  // State for axis selection and log scale toggles
  let axes = $state({
    x: ALL_METRICS.κ_SRME,
    y: ALL_METRICS.CPS,
    color_value: ALL_METRICS.F1,
    size_value: HYPERPARAMS.model_params,
  })
  let log = $state({ x: false, y: false, color_value: false, size_value: false })
  let is_fullscreen = $state(false)

  // State for extra controls
  let show_extra_controls = $state(false)
  let selected_color_schemes = $state<D3ColorSchemeName[]>([`Viridis`])
  let x_ticks = $state(5)
  let y_ticks = $state(5)
  let x_grid = $state(true)
  let y_grid = $state(true)

  const { model_params, graph_construction_radius, max_force, max_steps } = HYPERPARAMS
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
        // For each property (x, y, color, size), get the value based on the selected property path
        let x_path = `${axes.x?.path ?? ``}.${axes.x?.key ?? ``}`
        if (x_path.startsWith(`.`)) x_path = x_path.slice(1)
        let x_val = get_nested_value(model, x_path)
        if (x_path.includes(`date`)) x_val = new Date(x_val as string)

        let y_path = `${axes.y?.path ?? ``}.${axes.y?.key ?? ``}`
        if (y_path.startsWith(`.`)) y_path = y_path.slice(1)
        let y_val = get_nested_value(model, y_path)
        if (y_path.includes(`date`)) y_val = new Date(y_val as string)
        let color_path = `${axes.color_value?.path ?? ``}.${axes.color_value?.key ?? ``}`
        if (color_path.startsWith(`.`)) color_path = color_path.slice(1)
        let color_value = get_nested_value(model, color_path)
        if (color_path.includes(`date`) && color_value instanceof Date)
          color_value = new Date(color_value)

        let size_path = `${axes.size_value?.path ?? ``}.${axes.size_value?.key ?? ``}`
        if (size_path.startsWith(`.`)) size_path = size_path.slice(1)
        let size_value = get_nested_value(model, size_path)
        if (axes.size_value?.key === date_key) {
          const timestamp = new Date(String(size_value)).getTime()
          if (!isNaN(timestamp)) size_value = timestamp
        }

        // Prepare metadata for display in tooltip
        const { model_name, date_added, color } = model
        const metadata = {
          model_name,
          date_added,
          days_ago: calculate_days_ago(model.date_added),
        }
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

  // Update series when dependencies change
  let series: DataSeries[] = $derived.by(() => {
    // Base styling for points
    const point_styles = plot_data.map((item) => ({
      fill: point_color ?? item.color ?? `#4dabf7`,
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
      size_values: axes.size_value
        ? plot_data.map((item) => item.size_value as number)
        : undefined,
    }

    if (!show_model_labels) return [base_series]

    base_series.point_label = plot_data.map((item) => ({
      text: item.metadata.model_name,
      offset_y: 0,
      offset_x: 10,
      font_size: `14px`,
      color: `black`,
      auto_placement: true,
    }))

    return [base_series]
  })

  function handle_keydown(event: KeyboardEvent): void {
    if (event.key === `Escape`) {
      if (is_fullscreen) {
        is_fullscreen = false
        document.body.classList.remove(`fullscreen`)
      } else if (show_extra_controls) {
        show_extra_controls = false
      }
    }
  }
</script>

<svelte:window on:keydown={handle_keydown} />

<div class="plot-container" class:fullscreen={is_fullscreen}>
  <div class="top-buttons">
    <button
      class="fullscreen-toggle icon-button"
      onclick={() => {
        is_fullscreen = !is_fullscreen
        document.body.classList.toggle(`fullscreen`, is_fullscreen)
      }}
      aria-label={is_fullscreen ? `Exit fullscreen` : `Enter fullscreen`}
      title={is_fullscreen ? `Exit fullscreen` : `Enter fullscreen`}
    >
      <svg style="width: 1.3em; height: 1.3em;">
        <use href="#icon-{is_fullscreen ? `close` : `maximize`}" />
      </svg>
    </button>
    <button
      class="settings-toggle icon-button"
      onclick={() => (show_extra_controls = !show_extra_controls)}
      aria-label="Toggle extra plot controls"
      title="Toggle extra plot controls"
    >
      <svg style="width: 1.3em; height: 1.3em;">
        <use href="#icon-settings" />
      </svg>
    </button>
  </div>

  {#if show_extra_controls}
    <div
      use:click_outside={{ callback: () => (show_extra_controls = false) }}
      class="extra-controls"
    >
      <label style="grid-column: 1/-1">
        <input type="checkbox" bind:checked={show_model_labels} /> Show Labels
      </label>
      <label>
        <input type="checkbox" bind:checked={x_grid} /> X Grid
      </label>
      <label for="x-ticks">Ticks:</label>
      <input id="x-ticks" type="number" min="0" max="20" bind:value={x_ticks} />

      <label>
        <input type="checkbox" bind:checked={y_grid} /> Y Grid
      </label>
      <label for="y-ticks">Ticks:</label>
      <input id="y-ticks" type="number" min="0" max="20" bind:value={y_ticks} />
      <ColorScaleSelect
        bind:selected={selected_color_schemes}
        style="margin: 0; grid-column: 1/-1;"
      />
    </div>
  {/if}

  <div class="controls-grid">
    <!-- prettier-ignore -->
    {#each [
      { id: `x`, label: `X Axis`, log_state: log.x },
      { id: `y`, label: `Y Axis`, log_state: log.y },
      { id: `color_value`, label: `Color`, log_state: log.color_value },
      { id: `size_value`, label: `Size`, log_state: log.size_value },
    ] as const as control (control.id)}
      {@const [min, max] = extent(plot_data, (d) => d[control.id] as number)}
      <!-- want at least two orders of magnitude difference between min and max for log scaling to make sense -->
      {@const disable_log = Boolean(min === undefined || max === undefined || min <= 0 || 100 * min > max)}
      <label for={control.id}>{control.label}</label>
      <Select
        id={control.id}
        selected={[axes[control.id]]}
        bind:value={axes[control.id]}
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
        {@html format_property_path(`${prop.path ?? ``}.${prop.short ?? prop.label}`.replace(/^\./, ``))}
        <span style="font-size: smaller; color: gray; margin-left: 0.5em;">
          ({model_counts_by_prop[prop.key]} models)
        </span>
      </Select>
      <label
        aria-disabled={disable_log}
        style="transition: opacity 0.2s;"
        style:visibility={disable_log ? `hidden` : `visible`}
      >
        <input type="checkbox" bind:checked={log[control.id]} disabled={disable_log} />
        Log scale
      </label>
    {/each}
  </div>

  <div
    class:full-bleed-1400={!is_fullscreen}
    style="height: {is_fullscreen ? `100%` : `600px`}; margin-block: 1em;"
  >
    <!-- TODO fix x_lim and y_lim to use metric ranges-->
    <ScatterPlot
      {series}
      x_label="{axes.x?.svg_label ?? axes.x?.label} {axes.x?.better
        ? `<tspan style='font-size: 0.8em;'>(${axes.x?.better}=better)</tspan>`
        : ``}"
      y_label="{axes.y?.svg_label ?? axes.y?.label} {axes.y?.better
        ? `<tspan style='font-size: 0.8em;'>(${axes.y?.better}=better)</tspan>`
        : ``}"
      x_lim={axes.x?.range}
      y_lim={axes.y?.range}
      x_format={axes.x?.format}
      y_format={axes.y?.format}
      markers="points"
      {style}
      x_scale_type={log.x ? `log` : `linear`}
      y_scale_type={log.y ? `log` : `linear`}
      x_label_shift={{ y: -60 }}
      y_label_shift={{
        y: [date_key, params_key].includes(axes.y?.key ?? ``) ? -40 : -10,
      }}
      {x_ticks}
      {y_ticks}
      {x_grid}
      {y_grid}
      color_scale={{
        scheme: to_title(selected_color_schemes?.[0] ?? `Viridis`) as D3ColorSchemeName,
        type: log.color_value ? `log` : `linear`,
      }}
      size_scale={{ radius_range: [5, 20] }}
      color_bar={{
        title: `${axes.color_value?.label}${axes.color_value?.better ? ` (${axes.color_value?.better}=better)` : ``}`,
        margin: { t: 30, l: 80, b: 80, r: 50 },
        tick_format: axes.color_value?.format,
      }}
      label_placement_config={{
        link_strength: 5,
        link_distance: 5,
        link_distance_range: [15, 20],
      }}
      {...rest}
    >
      {#snippet tooltip({ x_formatted, y_formatted, metadata })}
        {#if metadata}
          {@const point = plot_data.find(
            (m) => m.metadata.model_name === metadata.model_name,
          )}
          <strong>{metadata.model_name}</strong><br />
          {@html axes.x.label}: {x_formatted}
          {#if axes.x.key === `date_added` && metadata.days_ago}
            <small>({metadata.days_ago} days ago)</small>{/if}<br />
          {@html axes.y.label}: {y_formatted}<br />
          {#if ![`model_params`, `date_added`].includes(axes.color_value?.key ?? ``) && point?.color_value !== undefined}
            {@html axes.color_value.label}:
            {pretty_num(point.color_value as number)}<br />
          {/if}
          {#if axes.size_value && point?.size_value !== undefined}
            {@html axes.size_value.label}:
            {pretty_num(point.size_value as number)}<br />
          {/if}
        {/if}
      {/snippet}
    </ScatterPlot>
  </div>
</div>

<style>
  /* Add styles globally to hide body scrollbar when a component is fullscreen */
  :global(body.fullscreen) {
    overflow: hidden;
  }
  .plot-container {
    position: relative; /* Needed for absolute positioning of the button */
    /* Use site background */
    background: var(--bg);
  }
  .plot-container.fullscreen {
    /* Ensure fixed positioning covers the entire viewport */
    position: fixed;
    top: 0;
    left: 0;
    width: 100vw;
    height: 100vh;
    z-index: 1000;
    /* Explicitly set background for fullscreen using the site's main background color */
    background: var(--night); /* Use --night from app.css */
    overflow: auto; /* Allow scrolling within the fullscreen container if content overflows */
    padding: 2em;
    box-sizing: border-box;
    /* Use flexbox for layout */
    display: flex;
    flex-direction: column;
    gap: 1em; /* Add gap between controls and plot */
  }
  .top-buttons {
    position: absolute;
    top: 0;
    right: -5em;
    display: flex;
    gap: 1ex;
  }
  .top-buttons button {
    cursor: pointer;
    padding: 8px;
    border-radius: 50%;
    transition: background-color 0.2s;
    background-color: rgba(255, 255, 255, 0.15);
  }
  .icon-button:hover {
    background-color: rgba(255, 255, 255, 0.3);
  }
  :global(body.fullscreen) .top-buttons {
    top: 1em;
    right: 1em;
  }
  .extra-controls {
    position: absolute;
    top: 45px; /* Adjust as needed to position below the gear icon */
    right: -5em; /* Align with the right edge of the gear icon */
    background-color: rgba(0, 0, 0, 0.8); /* Slightly transparent white */
    border-radius: 6px;
    padding: 5px 15px 15px;
    box-sizing: border-box;
    z-index: 2;
    display: grid;
    grid-template-columns: auto auto 1fr; /* Checkbox | Label | Input */
    gap: 5pt 1em;
  }
  :global(body.fullscreen) .extra-controls {
    top: 3.3em;
    right: 1em;
  }
  .extra-controls input[type='number'] {
    width: 30px;
  }
  .controls-grid {
    display: grid;
    grid-template-columns: auto 1fr auto; /* columns for: Label Select Checkbox */
    gap: 1ex;
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
