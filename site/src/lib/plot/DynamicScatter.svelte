<script module lang="ts">
  // Axis metadata is exported so pages can bind selections and render matching headings.
  import { ALL_METRICS, HYPERPARAMS, METADATA_COLS } from '$lib/labels'

  export const scatter_options = [
    ...Object.values(ALL_METRICS),
    HYPERPARAMS.model_params,
    METADATA_COLS.date_added,
    METADATA_COLS.n_training_materials,
    METADATA_COLS.n_training_structures,
    HYPERPARAMS.graph_construction_radius,
    HYPERPARAMS.max_force,
    HYPERPARAMS.max_steps,
    HYPERPARAMS.batch_size,
    HYPERPARAMS.epochs,
    HYPERPARAMS.n_layers,
  ]
  // Keyed lookup for bound axis/color selections.
  export const scatter_options_by_key = Object.fromEntries(
    scatter_options.map((opt) => [opt.key, opt]),
  )
  // Labels may contain HTML such as <sub>.
  export const scatter_axis_label = (key: string): string =>
    scatter_options_by_key[key]?.label ?? key
</script>

<script lang="ts">
  import { goto } from '$app/navigation'
  import type { ModelData } from '$lib'
  import { extent } from 'd3-array'
  import { format_num, format_relative_time, ScatterPlot } from 'matterviz'
  import type { D3InterpolateName } from 'matterviz/colors'
  import {
    ColorScaleSelect,
    create_color_scale,
    DEFAULT_SERIES_SYMBOLS,
  } from 'matterviz/plot'
  import type { DataSeries } from 'matterviz/plot'
  import type { ComponentProps } from 'svelte'
  import { tick } from 'svelte'
  import Select from 'svelte-multiselect'
  import { DISCOVERY_SET_LABELS, format_property_path } from '$lib/labels'
  import { get_nested_value, is_finite_num, label_data_path } from '$lib/metrics'
  import { make_models_legend } from '$lib/fig-helpers'
  import type { Label } from '$lib/types'

  // Keep size-select labels short by dropping discovery-set segments and abbreviating
  // "Geometry Optimization" to "Geo Opt".
  const discovery_set_keys = Object.keys(DISCOVERY_SET_LABELS)
  function format_size_option_path(path: string): string {
    const parts = path.split(`.`).filter((part) => !discovery_set_keys.includes(part))
    return format_property_path(parts.join(`.`)).replace(
      `Geometry Optimization`,
      `Geo Opt`,
    )
  }

  // Get value from model using label's path, converting dates to timestamps
  function get_label_value(model: ModelData, label: Label | undefined): unknown {
    const path = label_data_path(label)
    let val = get_nested_value(model, path)
    if (path.includes(`date`)) val = new Date(val as string).getTime()
    return val
  }

  const {
    legend_group,
    legend: models_legend,
    collapse_on_outside_click,
  } = make_models_legend()

  let {
    models,
    model_filter = () => true,
    point_color = null,
    show_model_labels = true,
    // Bindable so page headings can track user-selected axes.
    x_key = $bindable(ALL_METRICS.κ_SRME.key),
    y_key = $bindable(ALL_METRICS.CPS.key),
    color_key = $bindable(ALL_METRICS.F1.key),
    ...rest
  }: ComponentProps<typeof ScatterPlot> & {
    models: ModelData[]
    model_filter?: (model: ModelData) => boolean
    point_color?: string | null
    show_model_labels?: boolean
    x_key?: string
    y_key?: string
    color_key?: string
  } = $props()

  let log = $state({ x: false, y: false, color: false, size: false })
  let size_prop = $state(HYPERPARAMS.model_params as (typeof scatter_options)[number])

  let axes = $derived({
    x: scatter_options_by_key[x_key],
    y: scatter_options_by_key[y_key],
    color_value: scatter_options_by_key[color_key],
    size_value: size_prop,
  })

  let color_scheme: D3InterpolateName = $state(`interpolateViridis`)
  const ticks = 5
  let display = $state({ x_grid: true, y_grid: true })

  let size_multiplier = $state(1)
  let label_font_size = $state(12)
  let leader_line_threshold = $state(15)

  // Enable log scale only for positive ranges spanning at least two decades.
  const can_log = (ext: [number | undefined, number | undefined]): boolean =>
    ext[0] !== undefined && ext[0] > 0 && 100 * ext[0] <= (ext[1] ?? 0)

  let filtered_models = $derived(models.filter(model_filter))
  let model_counts_by_prop = $derived(
    Object.fromEntries(
      scatter_options.map((prop) => [
        prop.key,
        filtered_models.filter(
          (model) => get_nested_value(model, label_data_path(prop)) !== undefined,
        ).length,
      ]),
    ),
  )

  // Axis/color-select options with model counts.
  let prop_options = $derived(
    scatter_options.map((prop) => ({
      key: prop.key,
      label: `${prop.label} (${model_counts_by_prop[prop.key]} models)`,
      unit: prop.unit,
    })),
  )

  // ScatterPlot requests only x/y data here; the color bar handles color changes below.
  const data_loader = async (axis: string, key: string) => {
    if (axis === `x`) x_key = key
    else if (axis === `y`) y_key = key

    await tick()
    return { series, axis_label: scatter_axis_label(key) }
  }

  const format_label_title = (prop: Label | undefined): string =>
    `${prop?.label ?? ``}${prop?.better ? ` (${prop?.better}=better)` : ``}`
  const colorbar_tick_format = (prop: Label | undefined): string =>
    (prop?.format ?? `~s`).replace(/(?<precision>\.\d+)f$/, `$<precision>~f`)

  interface PointMetadata extends Record<string, unknown> {
    model_name: string
    date_added: string
    days_ago: string
    model_key?: string
  }

  let plot_data = $derived(
    filtered_models.flatMap((model) => {
      const x = get_label_value(model, axes.x)
      const y = get_label_value(model, axes.y)
      const color_value = get_label_value(model, axes.color_value)
      const size_value = get_label_value(model, axes.size_value)
      if (
        !is_finite_num(x) ||
        !is_finite_num(y) ||
        !is_finite_num(size_value) ||
        (point_color === null && !is_finite_num(color_value))
      ) {
        return []
      }

      const { model_name, date_added: model_date, model_key } = model
      const days_ago = format_relative_time(model.date_added)
      const metadata: PointMetadata = {
        model_name,
        date_added: model_date,
        days_ago,
        model_key,
      }
      return [{ x, y, color_value, size_value, metadata }]
    }),
  )

  // Mirror the plot color scale so built-in legend swatches match their points.
  let legend_color_scale = $derived.by(() => {
    const [min, max] = extent(plot_data, (item) => item.color_value as number)
    return create_color_scale(
      { scheme: color_scheme, type: log.color ? `log` : `linear` },
      [min ?? 0, max ?? 1],
    )
  })
  const point_fill = (color_value: unknown): string =>
    point_color ?? (legend_color_scale(color_value as number) as string) ?? `gray`

  let can_log_x = $derived(can_log(extent(plot_data, (d) => d.x)))
  let can_log_y = $derived(can_log(extent(plot_data, (d) => d.y)))
  let can_log_color = $derived(can_log(extent(plot_data, (d) => d.color_value as number)))
  let can_log_size = $derived(can_log(extent(plot_data, (d) => d.size_value)))

  // One series per model enables per-model legend toggles and distinct marker shapes.
  let series: DataSeries<PointMetadata>[] = $derived(
    plot_data.map((item, idx) => ({
      x: [item.x],
      y: [item.y],
      label: item.metadata.model_name,
      legend_group,
      markers: `points` as const,
      metadata: [item.metadata],
      point_style: {
        fill: point_fill(item.color_value),
        symbol_type: DEFAULT_SERIES_SYMBOLS[idx % DEFAULT_SERIES_SYMBOLS.length],
      },
      color_values: point_color === null ? [item.color_value as number] : undefined,
      size_values: axes.size_value ? [item.size_value] : undefined,
      point_label: show_model_labels
        ? [
            {
              text: item.metadata.model_name,
              font_size: `${label_font_size}px`,
              auto_placement: true,
            },
          ]
        : [],
    })),
  )
</script>

<div
  class="bleed-1400 collapsible-legend"
  style="margin-block: 2em"
  {@attach collapse_on_outside_click}
>
  <div class="controls-row">
    <label for="size-select">Marker Size</label>
    <Select
      options={scatter_options}
      id="size-select"
      bind:value={size_prop}
      maxSelect={1}
      minSelect={1}
      style="flex: 1; max-width: 300px; margin: 0"
      ulSelectedStyle="flex-wrap: nowrap; overflow: hidden; min-width: 0;"
      liSelectedStyle="font-size: 14px; min-width: 0; max-width: 100%; overflow: hidden;"
      liOptionStyle="font-size: 13px;"
    >
      {#snippet children({
        option: prop,
        type,
      }: {
        option: (typeof scatter_options)[number]
        type: string
      })}
        <span class:selected-label={type === `selected`}>
          {@html format_size_option_path(label_data_path(prop))}
          <span style="font-size: smaller; color: gray">
            {model_counts_by_prop[prop.key]} models
          </span>
        </span>
      {/snippet}
    </Select>
  </div>

  <ScatterPlot
    style="height: 600px"
    bind:series
    legend={models_legend}
    padding={{ b: 70 }}
    x_axis={{
      label: axes.x?.label,
      format: axes.x?.format,
      scale_type: log.x ? `log` : `linear`,
      label_shift: { y: -50 },
      ticks,
      options: prop_options,
      selected_key: x_key,
    }}
    y_axis={{
      label: axes.y?.label,
      format: axes.y?.format,
      scale_type: log.y ? `log` : `linear`,
      label_shift: {
        x: -10,
        y: [`date_added`, `model_params`].includes(axes.y?.key ?? ``) ? -40 : -10,
      },
      ticks,
      options: prop_options,
      selected_key: y_key,
    }}
    bind:display
    color_scale={{ scheme: color_scheme, type: log.color ? `log` : `linear` }}
    size_scale={{
      radius_range: [5 * size_multiplier, 10 * size_multiplier],
      type: log.size ? `log` : `linear`,
    }}
    color_bar={{
      title: format_label_title(axes.color_value),
      margin: { t: 30, l: 80, b: 80, r: 50 },
      tick_format: colorbar_tick_format(axes.color_value),
      property_options: prop_options,
      selected_property_key: color_key,
      data_loader: async (key) => {
        color_key = key
        const prop = scatter_options_by_key[key]
        const values = filtered_models
          .map((model) => get_label_value(model, prop))
          .filter((val): val is number => typeof val === `number` && isFinite(val))
        const [min, max] = extent(values)
        return {
          range: [min ?? 0, max ?? 1],
          title: format_label_title(prop),
        }
      },
    }}
    label_placement_config={{
      leader_line_threshold,
      max_neighbors: { count: 3, radius: 40 },
    }}
    point_events={{
      onclick: ({ point }) => goto(`/models/${point.metadata?.model_key ?? ``}`),
    }}
    {...rest}
    {data_loader}
  >
    {#snippet controls_extra()}
      <div class="log-toggles" style="display: flex; gap: 1em; flex-wrap: wrap">
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

      <label title="Toggle visibility of model name labels on the scatter plot points">
        <input type="checkbox" bind:checked={show_model_labels} /> Show Labels
      </label>
      <ColorScaleSelect bind:value={color_scheme} style="margin: 0" />
      <label
        for="size-multiplier"
        title="Adjust the base size of all points on the scatter plot (multiplier for radius)"
        >Point Size
      </label>
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
        >Label Size
      </label>
      <input
        id="label-font-size"
        type="range"
        min="8"
        max="24"
        step="1"
        bind:value={label_font_size}
      />
      <label
        title="Minimum label displacement in pixels before drawing a leader line"
        for="leader-line-threshold">Leader Line</label
      >
      <div class="combined-link-controls">
        <input
          id="leader-line-threshold"
          type="number"
          min="0"
          max="100"
          bind:value={leader_line_threshold}
          title="Leader line threshold"
        />
      </div>
    {/snippet}

    {#snippet tooltip({ x_formatted, y_formatted, metadata })}
      {#if metadata}
        {@const point = plot_data.find(
          (m) => m.metadata.model_name === metadata.model_name,
        )}
        <strong>{metadata.model_name}</strong><br />
        {@html axes.x?.label}: {x_formatted}
        {#if axes.x?.key === `date_added` && metadata.days_ago}
          <small>({metadata.days_ago})</small>{/if}<br />
        {@html axes.y?.label}: {y_formatted}<br />
        {#if ![`model_params`, `date_added`].includes(axes.color_value?.key ?? ``) && point?.color_value !== undefined}
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
    gap: 1ex 0.6em;
    /* asymmetric margin counterbalances the collapsed legend parked right of center
       (a `translate` would be simpler but turns the row into a containing block for
       the marker-size dropdown's position: fixed options list, breaking it) */
    margin: 0 12em 1em 3em;
    justify-content: center;
    /* paint above the (later-DOM) plot so the open dropdown stays interactive */
    position: relative;
    z-index: 1;
  }
  div.controls-row label {
    font-weight: 500;
    font-size: 14px;
  }
  /* move ScatterPlot's legend into the controls row */
  div.bleed-1400 :global(.scatter > .legend) {
    top: -42px !important;
    bottom: auto !important;
    font-size: 14px;
  }
  /* collapsed: show only the group header beside the size select */
  div.bleed-1400 :global(.scatter > .legend:not(:has(.legend-item))) {
    left: calc(50% + 10em) !important;
    width: max-content !important;
    justify-content: flex-start !important;
  }
  /* expanded: wrap model items across the plot width */
  div.bleed-1400 :global(.scatter > .legend:has(.legend-item)) {
    left: 10px !important;
    width: calc(100% - 20px) !important;
  }
  div.combined-link-controls {
    display: flex;
    gap: 0.5em;
    align-items: center;
  }
  div.combined-link-controls input[type='number'] {
    width: 50px;
  }
  span.selected-label {
    display: block;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
  }
</style>
