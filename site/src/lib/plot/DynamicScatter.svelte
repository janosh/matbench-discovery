<script lang="ts" generics="T extends object">
  import { goto } from '$app/navigation'
  import type { Label, DiscoverySet } from '$lib/types'
  import { extent } from 'd3-array'
  import { format_value_or_num } from 'matterviz/labels'
  import { ScatterPlot } from 'matterviz/plot'
  import type { AxisConfig, DataSeries, InternalPoint } from 'matterviz/plot'
  import type { ComponentProps } from 'svelte'
  import { MultiSelect } from 'svelte-widgets'
  import {
    ALL_METRICS,
    DISCOVERY_SET_LABELS,
    format_property_path,
    format_relative_time,
    HYPERPARAMS,
    scatter_options,
    title_case,
  } from '$lib/labels'
  import {
    get_nested_value,
    is_finite_num,
    label_data_path,
    metric_data_path,
  } from '$lib/metrics'
  import { make_models_legend } from '$lib/fig-helpers'
  import { pareto_staircase, sota_frontier_indices, sota_step_line } from '$lib/sota'

  // Keep size-select labels short by dropping discovery-set segments and abbreviating
  // "Geometry Optimization" to "Geo Opt".
  const discovery_set_keys = Object.keys(DISCOVERY_SET_LABELS)
  const format_size_option_path = (path: string): string =>
    format_property_path(
      path
        .split(`.`)
        .filter((part) => !discovery_set_keys.includes(part))
        .join(`.`),
    ).replace(`Geometry Optimization`, `Geo Opt`)

  type ScatterOption = Label &
    Record<string, unknown> & { categories?: Record<string, string> }

  // Resolve each axis path/date conversion once per selection, not per model.
  const label_accessor = (label: ScatterOption | undefined) => {
    const path = label ? metric_data_path(label, discovery_set) : ``
    const is_date = path.includes(`date`)
    const categories = label?.categories && Object.keys(label.categories)
    return (model: T): unknown => {
      const value = get_nested_value(model, path)
      if (value == null) return undefined
      if (categories) {
        const idx = categories.indexOf(String(value))
        return idx === -1 ? undefined : idx
      }
      return is_date ? new Date(String(value)).getTime() : value
    }
  }

  const {
    legend_group,
    legend: models_legend,
    collapse_on_outside_click,
    toggle: toggle_models,
  } = make_models_legend()

  let {
    models,
    model_filter = () => true,
    item_name = `models`,
    get_identity = (model: T) => ({
      key: String(get_nested_value(model, `model_key`)),
      name: String(get_nested_value(model, `model_name`)),
      href: `/models/${get_nested_value(model, `model_key`)}`,
    }),
    show_model_labels = true,
    // Bindable so page headings can track user-selected axes.
    x_key = $bindable(ALL_METRICS.κ_SRME.key),
    y_key = $bindable(ALL_METRICS.CPS.key),
    color_key = $bindable(ALL_METRICS.F1.key),
    options = scatter_options,
    show_pareto_frontier = false,
    highlight_keys,
    discovery_set = `unique_prototypes`,
    size_key = $bindable(HYPERPARAMS.model_params.key),
    legend = models_legend,
    bleed = true,
    ...rest
  }: ComponentProps<typeof ScatterPlot> & {
    models: T[]
    model_filter?: (model: T) => boolean
    item_name?: string
    get_identity?: (record: T) => { key: string; name: string; href: string }
    show_model_labels?: boolean
    x_key?: string
    y_key?: string
    color_key?: string
    // Labels selectable for the axes, color and size; keys must be unique
    options?: ScatterOption[]
    // trace the staircase of non-dominated models (needs better-direction on both axes)
    show_pareto_frontier?: boolean
    // when given, only these models are labeled and drawn at full opacity with a ring,
    // the rest recede into a translucent field
    highlight_keys?: Set<string>
    discovery_set?: DiscoverySet
    size_key?: string
    // span the full viewport width (the default on task pages, off inside dialogs)
    bleed?: boolean
  } = $props()

  const log_dims = [`x`, `y`, `color`, `size`] as const
  const log_dim_labels = { x: `X`, y: `Y`, color: `Color`, size: `Size` } as const
  let options_by_key = $derived(Object.fromEntries(options.map((opt) => [opt.key, opt])))
  let axes = $derived({
    x: options_by_key[x_key],
    y: options_by_key[y_key],
    color_value: options_by_key[color_key],
    size_value: options_by_key[size_key],
  })
  const color_entries = $derived(
    axes.color_value?.categories && Object.entries(axes.color_value.categories),
  )

  let axis_accessors = $derived(Object.values(axes).map(label_accessor))

  let display = $state({ x_grid: true, y_grid: true })

  let filtered_models = $derived(models.filter(model_filter))
  let models_by_name = $derived(
    Object.groupBy(filtered_models, (model) => get_identity(model).name),
  )
  const model_count = (prop: ScatterOption): number => {
    const accessor = label_accessor(prop)
    return filtered_models.filter((model) => is_finite_num(accessor(model))).length
  }

  const format_label_title = (prop: Label | undefined): string =>
    `${prop?.label ?? ``}${prop?.better ? ` (${prop?.better}=better)` : ``}`
  // fallback is a trimmed float, not `~s`: SI prefixes render 0.5 as "500m", and all
  // labels with big-count values (model params, training size) set format `~s` anyway
  const colorbar_tick_format = (prop: Label | undefined): string =>
    (prop?.format ?? `.2~f`).replace(/(?<precision>\.\d+)f$/, `$<precision>~f`)

  let plot_data = $derived(
    filtered_models.flatMap((model) => {
      const values = axis_accessors.map((accessor) => accessor(model))
      if (!values.every(is_finite_num)) return []
      const [x, y, color_value, size_value] = values
      const { name: model_name, key: model_key, href } = get_identity(model)
      const benchmark_added = get_nested_value(model, `dates.benchmark_added`)
      const days_ago =
        typeof benchmark_added === `string` ? format_relative_time(benchmark_added) : ``
      const metadata = {
        ...model,
        model_key,
        model_name,
        href,
        days_ago,
        color_value,
        size_value,
      }
      return [{ x, y, color_value, size_value, metadata }]
    }),
  )
  type PointMetadata = (typeof plot_data)[number][`metadata`]

  // Log scales need positive values spanning at least two decades.
  const supports_log = (
    prop: Label | undefined,
    value_key: `x` | `y` | `color_value` | `size_value`,
  ): boolean => {
    const [min, max] = extent(plot_data, (point) => point[value_key])
    return (
      !label_data_path(prop).includes(`date`) &&
      min !== undefined &&
      min > 0 &&
      100 * min <= max
    )
  }
  let can_log = $derived({
    x: supports_log(axes.x, `x`),
    y: supports_log(axes.y, `y`),
    color: supports_log(axes.color_value, `color_value`),
    size: supports_log(axes.size_value, `size_value`),
  })
  let supported_log_dims = $derived(log_dims.filter((dim) => can_log[dim]))
  // Initialize automatically; manual toggles reset when data or dimensions change.
  let log = $derived({ ...can_log })
  const scale_of = (dim: keyof typeof log) =>
    log[dim] ? (`log` as const) : (`linear` as const)

  // Staircase through the non-dominated models, tracing the boundary of the dominated
  // region. With a date on the x-axis it becomes the running best over time (records
  // extended to today); otherwise it needs a better-direction on both axes.
  let pareto_series = $derived.by((): DataSeries<PointMetadata> | null => {
    const [x_better, y_better] = [axes.x?.better, axes.y?.better]
    if (!show_pareto_frontier || !y_better) return null

    const line_series = (label: string, xs: number[], ys: number[]) => ({
      x: xs,
      y: ys,
      label,
      legend_group,
      markers: `line` as const,
      line_style: { stroke: `gray`, stroke_width: 1.5, line_dash: `5 3` },
    })

    if (label_data_path(axes.x).includes(`date`)) {
      // field-progress view: which releases moved the frontier, and where it stands
      const points = plot_data.map((pt) => ({ date: pt.x, value: pt.y }))
      const records = sota_frontier_indices(points, y_better).map((idx) => points[idx])
      if (records.length === 0) return null
      const { x, y } = sota_step_line(records, Date.now())
      return line_series(`Running best`, x, y)
    }

    if (!x_better) return null
    const staircase = pareto_staircase(plot_data, x_better, y_better)
    return staircase && line_series(`Pareto frontier`, staircase.x, staircase.y)
  })

  // Suppress hover on the frontier line: its staircase corners are not models (no
  // metadata), so snapping to them showed a useless tooltip. Frontier vertices that ARE
  // models still get the model tooltip: the model series win the closest-point tie by
  // coming first in `series`.
  // Axis or data changes invalidate the previous point's values and category indices.
  let tooltip_point: InternalPoint | null = $derived.by(() => {
    void plot_data
    return null
  })
  $effect(() => {
    if (tooltip_point && !tooltip_point.metadata) tooltip_point = null
  })

  // without highlight_keys every model is in focus
  const is_dimmed = ({ model_key }: PointMetadata) =>
    highlight_keys !== undefined && !highlight_keys.has(model_key)
  // One series per model enables per-model legend toggles. Highlighted models come last
  // so they paint over the dimmed field.
  let series: DataSeries<PointMetadata>[] = $derived([
    ...plot_data
      .toSorted(
        (pt1, pt2) => Number(is_dimmed(pt2.metadata)) - Number(is_dimmed(pt1.metadata)),
      )
      .map(({ x, y, color_value, size_value, metadata }) => {
        const { model_key, model_name } = metadata
        const dimmed = is_dimmed(metadata)
        return {
          id: model_key,
          x: [x],
          y: [y],
          label:
            (models_by_name[model_name]?.length ?? 0) > 1
              ? `${model_name} (${model_key})`
              : model_name,
          legend_group,
          markers: `points` as const,
          metadata,
          // uniform circles: color and size already encode data, and cycling 7 shapes
          // across 30+ models distinguished nothing while adding visual noise
          point_style: {
            fill: color_entries?.[color_value]?.[1],
            symbol_type: `Circle` as const,
            ...(dimmed
              ? { fill_opacity: 0.3 }
              : highlight_keys && { stroke: `currentColor`, stroke_width: 1.5 }),
          },
          color_values: color_entries ? undefined : [color_value],
          size_values: [size_value],
          point_label:
            show_model_labels && !dimmed
              ? { text: model_name, font_size: `12px`, auto_placement: true }
              : undefined,
        }
      }),
    ...(pareto_series ? [pareto_series] : []),
  ])

  const numeric_options = $derived(options.filter((option) => !option.categories))
  const axis_config = (dim: 'x' | 'y'): AxisConfig => ({
    options: numeric_options,
    selected_key: axes[dim]?.key,
    label: axes[dim]?.label,
    format: axes[dim]?.format,
    scale_type: label_data_path(axes[dim]).includes(`date`) ? `time` : scale_of(dim),
    ticks: 5,
  })
  const picker_id = $props.id()
</script>

<div
  class={[`dynamic-scatter collapsible-legend`, bleed && `bleed-1400`]}
  style="margin-block: 2em"
  {@attach collapse_on_outside_click}
>
  <div class="controls-row">
    <div class="property-picker">
      <label for="{picker_id}-size">Marker size</label>
      <MultiSelect
        options={numeric_options}
        id={`${picker_id}-size`}
        value={axes.size_value}
        mode="single"
        min_select={1}
        key={(opt: ScatterOption) => opt.key}
        on_change={(event) => {
          if (event.type === `add` && event.option) size_key = event.option.key
        }}
        style="flex: 1; min-width: 0; margin: 0; --sms-min-height: 28px"
        ul_selected_style="flex-wrap: nowrap; overflow: hidden; min-width: 0;"
        li_selected_style="font-size: 14px; min-width: 0; max-width: 100%; overflow: hidden;"
      >
        {#snippet children({ option: prop, type })}
          <span class:selected-label={type === `selected`}>
            {@html prop.label}
            <small
              >{format_size_option_path(label_data_path(prop))} · {model_count(prop)}
              {item_name}</small
            >
          </span>
        {/snippet}
      </MultiSelect>
    </div>
    {#if supported_log_dims.length}
      <div class="log-controls" role="group" aria-label="Logarithmic scales">
        <strong>Log Scale</strong>
        {#each supported_log_dims as dim (dim)}
          <label>
            <input
              type="checkbox"
              checked={log[dim]}
              onchange={(event) => (log = { ...log, [dim]: event.currentTarget.checked })}
            />
            {log_dim_labels[dim]}
          </label>
        {/each}
      </div>
    {/if}
    {#if legend && models_legend.collapsed_groups?.has(legend_group)}
      <button
        type="button"
        class="models-toggle"
        aria-expanded="false"
        onclick={toggle_models}
      >
        ▶ {legend_group}
      </button>
    {/if}
  </div>

  {#if !plot_data.length}
    <p style="text-align: center">No {item_name} have values for all selected columns.</p>
  {/if}

  <ScatterPlot
    style="height: 600px"
    {series}
    bind:tooltip_point
    {legend}
    padding={{ b: 70 }}
    x_axis={axis_config(`x`)}
    y_axis={axis_config(`y`)}
    on_axis_change={(axis, key) => {
      if (axis === `x`) x_key = key
      else if (axis === `y`) y_key = key
    }}
    bind:display
    color_scale={{ scheme: `interpolateViridis`, type: scale_of(`color`) }}
    size_scale={{
      radius_range: [5, 10],
      type: scale_of(`size`),
    }}
    color_bar={{
      title: format_label_title(axes.color_value),
      tick_format: colorbar_tick_format(axes.color_value),
      property_options: options,
      selected_property_key: color_key,
      on_property_change: (key) => (color_key = key),
      categories:
        color_entries &&
        Object.fromEntries(
          color_entries.map(([category, color]) => [title_case(category), color]),
        ),
    }}
    label_placement_config={{
      leader_line_threshold: 15,
      max_neighbors: { count: 3, radius: 40 },
      // Avoid blocking point-tween frames while axes change scale.
      sa_iterations: 100,
    }}
    point_events={{
      onclick: ({ point }) => {
        const href = point.metadata?.href
        if (typeof href === `string`) goto(href)
      },
    }}
    {...rest}
  >
    {#snippet controls_extra()}
      <label title="Toggle point labels on the plot">
        <input type="checkbox" bind:checked={show_model_labels} /> Show Labels
      </label>
    {/snippet}

    {#snippet tooltip({ x_formatted, y_formatted, metadata })}
      {#if metadata}
        <strong>{metadata.model_name}</strong><br />
        {@html axes.x?.label}: {x_formatted}
        {#if axes.x?.key === `benchmark_added` && metadata.days_ago}
          <small>({metadata.days_ago})</small>{/if}<br />
        {@html axes.y?.label}: {y_formatted}<br />
        {#if is_finite_num(metadata.color_value) && ![`model_params`, `benchmark_added`].includes(axes.color_value?.key ?? ``)}
          {@html axes.color_value?.label}:
          {color_entries
            ? title_case(color_entries[metadata.color_value][0])
            : format_value_or_num(metadata.color_value, axes.color_value?.format)}<br />
        {/if}
        {#if is_finite_num(metadata.size_value)}
          {@html axes.size_value.label}:
          {format_value_or_num(metadata.size_value, axes.size_value.format)}<br />
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
    justify-content: center;
    gap: 1ex 0.6em;
    margin: 0 0 1em;
    /* paint above the (later-DOM) plot so the open dropdown stays interactive */
    position: relative;
    z-index: 1;
  }
  div.log-controls,
  div.log-controls label,
  .property-picker {
    display: flex;
    align-items: center;
  }
  div.log-controls {
    gap: 0.6em;
    font-size: 14px;
    white-space: nowrap;
  }
  div.log-controls label {
    gap: 0.25em;
  }
  div.controls-row label {
    font-weight: 500;
    font-size: 14px;
  }
  button.models-toggle {
    padding: 2px 4px;
    background: none;
    color: inherit;
    font: inherit;
    font-size: 14px;
    white-space: nowrap;
  }
  /* align ScatterPlot's expanded legend with the controls row */
  div.dynamic-scatter :global(.scatter > .legend) {
    top: -42px !important;
    bottom: auto !important;
    font-size: 14px;
  }
  /* the controls-row button replaces the collapsed legend shell */
  div.dynamic-scatter:has(button.models-toggle) :global(.scatter > .legend) {
    display: none !important;
  }
  /* expanded: wrap model items across the plot width */
  div.dynamic-scatter :global(.scatter > .legend:has(.legend-item)) {
    left: 10px !important;
    width: calc(100% - 20px) !important;
  }
  span.selected-label {
    display: block;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
  }
  .property-picker {
    gap: 0.5em;
    flex: 1 1 220px;
    min-width: 0;
    max-width: 360px;
    label {
      white-space: nowrap;
    }
  }
  .selected-label small {
    display: none;
  }
</style>
