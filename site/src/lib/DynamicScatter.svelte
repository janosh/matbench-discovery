<script lang="ts">
  import { goto } from '$app/navigation'
  import type { ModelData } from '$lib'
  import { calculate_days_ago, Icon } from '$lib'
  import { extent } from 'd3-array'
  import { ColorScaleSelect, format_num, ScatterPlot } from 'matterviz'
  import Select from 'svelte-multiselect'
  import { click_outside, titles_as_tooltips } from 'svelte-zoo'
  import {
    ALL_METRICS,
    format_property_path,
    HYPERPARAMS,
    METADATA_COLS,
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

  // State for plot controls
  let show_controls = $state(false)
  let color_scheme = $state(`Viridis`)
  let x_ticks = $state(5)
  let y_ticks = $state(5)
  let x_grid = $state(true)
  let y_grid = $state(true)

  // State for point size and label placement
  let size_multiplier = $state(1)
  let label_font_size = $state(14) // Default font size
  let link_strength = $state(5)
  let link_distance = $state(5)
  let min_link_distance = $state(15)
  let max_link_distance = $state(20)

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
  ] as const
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
    if (typeof val === `number` && !isNaN(val)) return true
    if (val instanceof Date) return !isNaN(val.getTime()) // Check for valid dates
    return false
  }
  let plot_data = $derived(
    models
      .filter(model_filter)
      .map((model) => {
        // For each property (x, y, color, size), get the value based on the selected property path
        let x_path = `${axes.x?.path ?? ``}.${axes.x?.key ?? ``}`
        if (x_path.startsWith(`.`)) x_path = x_path.slice(1)
        let x_val = get_nested_value(model, x_path)
        if (x_path.includes(`date`)) x_val = new Date(x_val as string).getTime()

        let y_path = `${axes.y?.path ?? ``}.${axes.y?.key ?? ``}`
        if (y_path.startsWith(`.`)) y_path = y_path.slice(1)
        let y_val = get_nested_value(model, y_path)
        if (y_path.includes(`date`)) y_val = new Date(y_val as string).getTime()
        let color_path = `${axes.color_value?.path ?? ``}.${
          axes.color_value?.key ?? ``
        }`
        if (color_path.startsWith(`.`)) color_path = color_path.slice(1)
        let color_value = get_nested_value(model, color_path)
        if (color_path.includes(`date`)) {
          color_value = new Date(color_value as string).getTime()
        }

        let size_path = `${axes.size_value?.path ?? ``}.${axes.size_value?.key ?? ``}`
        if (size_path.startsWith(`.`)) size_path = size_path.slice(1)
        let size_value = get_nested_value(model, size_path)
        if (axes.size_value?.key === date_key) {
          const timestamp = new Date(String(size_value)).getTime()
          if (!isNaN(timestamp)) size_value = timestamp
        }

        // Prepare metadata for display in tooltip
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
    point_style: plot_data.map((item) => ({ fill: point_color ?? item.color })),
    metadata: plot_data.map((item) => item.metadata),
    color_values: point_color === null
      ? (plot_data
        .map((item) => item.color_value) // Map to potential Date or number
        .filter((val): val is number | Date => val !== undefined) // Filter out undefined
        .map((val) => (val instanceof Date ? val.getTime() : val)) as number[]) // Convert Dates to timestamps
      : undefined,
    size_values: axes.size_value
      ? plot_data.map((item) => item.size_value as number)
      : undefined,
    point_label: (show_model_labels ? plot_data : []).map((item) => ({
      text: item.metadata.model_name,
      offset_y: 0,
      offset_x: 10,
      font_size: `${label_font_size}px`, // Use state variable
      color: `black`,
      auto_placement: true,
    })),
  })

  function handle_keydown(event: KeyboardEvent): void {
    if (event.key === `Escape`) {
      if (is_fullscreen) {
        is_fullscreen = false
        document.body.classList.remove(`fullscreen`)
      } else show_controls = false
    }
  }

  interface DraggableOptions {
    handle_selector?: string
  }
  function draggable(node: HTMLElement, options?: DraggableOptions) {
    let dragging = false
    let start = { x: 0, y: 0 }
    let initial = { left: 0, top: 0, width: 0 }

    const handle = options?.handle_selector
      ? node.querySelector<HTMLElement>(options.handle_selector)
      : node

    if (!handle) return // Handle not found

    function handle_mousedown(event: MouseEvent) {
      // Only drag if mousedown is on the handle itself
      if (event.target !== handle) return

      dragging = true
      initial.left = node.offsetLeft
      initial.top = node.offsetTop
      initial.width = node.offsetWidth
      node.style.left = `${initial.left}px`
      node.style.top = `${initial.top}px`
      node.style.width = `${initial.width}px`
      node.style.right = `auto` // Prevent conflict with left
      start = { x: event.clientX, y: event.clientY }
      document.body.style.userSelect = `none` // Prevent text selection during drag
      handle!.style.cursor = `grabbing`
      window.addEventListener(`mousemove`, handle_mousemove)
      window.addEventListener(`mouseup`, handle_mouseup)
    }

    function handle_mousemove(event: MouseEvent) {
      if (!dragging) return
      const dx = event.clientX - start.x
      const dy = event.clientY - start.y
      node.style.left = `${initial.left + dx}px`
      node.style.top = `${initial.top + dy}px`
    }

    function handle_mouseup(event: MouseEvent) {
      dragging = false
      event.stopPropagation()
      document.body.style.userSelect = ``
      handle!.style.cursor = `grab`
      window.removeEventListener(`mousemove`, handle_mousemove)
      window.removeEventListener(`mouseup`, handle_mouseup)
    }

    handle.addEventListener(`mousedown`, handle_mousedown)
    handle.style.cursor = `grab` // Set initial cursor on handle

    return {
      destroy() {
        handle.removeEventListener(`mousedown`, handle_mousedown)
        window.removeEventListener(`mousemove`, handle_mousemove) // Clean up just in case
        window.removeEventListener(`mouseup`, handle_mouseup)
        handle.style.cursor = `` // Reset cursor
      },
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
      <Icon
        icon={is_fullscreen ? `Close` : `Maximize`}
        style="width: 1.3em; height: 1.3em"
      />
    </button>
    <button
      class="settings-toggle icon-button"
      onclick={() => (show_controls = !show_controls)}
      aria-label="Toggle plot controls"
    >
      <Icon icon="Settings" style="width: 1.3em; height: 1.3em" />
    </button>
  </div>

  {#if show_controls}
    <div
      use:click_outside={{ callback: () => (show_controls = false) }}
      use:draggable={{ handle_selector: `.drag-handle` }}
      class="controls"
      use:titles_as_tooltips
    >
      <Icon
        icon="DragIndicator"
        class="drag-handle"
        style="width: 1.3em; height: 1.3em; position: absolute; top: 5px; right: 5px; background-color: rgba(255, 255, 255, 0.2); border-radius: 3px"
      />
      <label
        style="grid-column: 1/-1"
        title="Toggle visibility of model name labels on the scatter plot points"
      >
        <input type="checkbox" bind:checked={show_model_labels} /> Show Labels
      </label>
      <ColorScaleSelect
        bind:value={color_scheme}
        selected={[color_scheme]}
        style="margin: 0; grid-column: 1/-1"
      />
      <label
        style="grid-column: 1 / -1"
        title="Toggle the visibility of vertical grid lines and set the approximate number of ticks on the X axis"
      >
        <input type="checkbox" bind:checked={x_grid} /> X Grid
        <label for="x-ticks" style="margin-left: 1em">Ticks:</label>
        <input
          id="x-ticks"
          type="number"
          min="0"
          max="20"
          bind:value={x_ticks}
          style="width: 50px"
        />
      </label>

      <label
        style="grid-column: 1 / -1"
        title="Toggle the visibility of horizontal grid lines and set the approximate number of ticks on the Y axis"
      >
        <!-- Span both columns -->
        <input type="checkbox" bind:checked={y_grid} /> Y Grid
        <label for="y-ticks" style="margin-left: 1em">Ticks:</label>
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
          style="flex-: 1"
        />
      </div>
    </div>
  {/if}

  <div class="controls-grid">
    {#each [
        { id: `x`, label: `X Axis`, log_state: log.x },
        { id: `y`, label: `Y Axis`, log_state: log.y },
        { id: `color_value`, label: `Color`, log_state: log.color_value },
        { id: `size_value`, label: `Size`, log_state: log.size_value },
      ] as const as
      control
      (control.id)
    }
      {@const [min, max] = extent(plot_data, (d) => d[control.id] as number)}
      <!-- want at least two orders of magnitude difference between min and max for log scaling to make sense -->
      {@const disable_log = Boolean(
        min === undefined || max === undefined || min <= 0 || 100 * min > max,
      )}
      <label for={control.id}>{control.label}</label>
      <Select
        {options}
        id={control.id}
        selected={[axes[control.id]]}
        bind:value={axes[control.id]}
        placeholder="Select {control.label}"
        maxSelect={1}
        minSelect={1}
        style="width: 100%; max-width: none; margin: 0"
        liSelectedStyle="font-size: 16px;"
        ulSelectedStyle="padding: 0;"
        --sms-selected-bg="none"
        --sms-border="1px solid rgba(255, 255, 255, 0.15)"
      >
        {#snippet children({ option: prop })}
          {@html format_property_path(
          `${prop.path ?? ``}.${prop.label ?? prop.label}`.replace(/^\./, ``),
        )}
          <span style="font-size: smaller; color: gray; margin-left: 0.5em">
            ({model_counts_by_prop[prop.key]} models)
          </span>
        {/snippet}
      </Select>
      <label
        aria-disabled={disable_log}
        style="transition: opacity 0.2s"
        style:visibility={disable_log ? `hidden` : `visible`}
      >
        <input type="checkbox" bind:checked={log[control.id]} disabled={disable_log} />
        Log scale
      </label>
    {/each}
  </div>

  <div
    class:full-bleed-1400={!is_fullscreen}
    style="margin-block: 1em"
    style:height={is_fullscreen ? `100%` : `600px`}
  >
    <ScatterPlot
      series={[series]}
      x_label="{axes.x?.svg_label ?? axes.x?.label} {axes.x?.better
        ? `<tspan style='font-size: 0.8em;'>(${axes.x?.better}=better)</tspan>`
        : ``}"
      y_label="{axes.y?.svg_label ?? axes.y?.label} {axes.y?.better
        ? `<tspan style='font-size: 0.8em;'>(${axes.y?.better}=better)</tspan>`
        : ``}"
      x_lim={axes.x?.range ? [...axes.x.range] : undefined}
      y_lim={axes.y?.range ? [...axes.y.range] : undefined}
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
      color_scale={{ scheme: color_scheme, type: log.color_value ? `log` : `linear` }}
      size_scale={{
        radius_range: [5 * size_multiplier, 20 * size_multiplier],
        type: log.size_value ? `log` : `linear`,
      }}
      color_bar={{
        title: `${axes.color_value?.label}${
          axes.color_value?.better ? ` (${axes.color_value?.better}=better)` : ``
        }`,
        margin: { t: 30, l: 80, b: 80, r: 50 },
        tick_format: axes.color_value?.format,
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
          {#if ![`model_params`, `date_added`].includes(axes.color_value?.key ?? ``) &&
          point?.color_value !== undefined}
            {@html axes.color_value.label}:
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
  .controls {
    position: absolute;
    top: 45px; /* Adjust as needed to position below the gear icon */
    right: -5em;
    max-width: 450px; /* Allow content width up to a max */
    background-color: var(--night, #1a1a1a);
    border: 1px solid rgba(255, 255, 255, 0.15);
    border-radius: 6px;
    padding: 5px 15px 15px;
    box-sizing: border-box;
    z-index: 2;
    display: grid;
    grid-template-columns: auto 1fr; /* Label | Control Area */
    gap: 8pt 1em;
  }
  .controls > *:nth-child(odd):not(label) {
    grid-column: 2;
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
  .combined-link-controls {
    display: flex;
    align-items: center;
    gap: 0.3em;
  }
  :global(body.fullscreen) .controls {
    top: 3.3em;
    right: 1em;
    left: auto; /* Ensure left is not set */
  }
  .controls input[type='number'] {
    width: 40px;
  }
</style>
