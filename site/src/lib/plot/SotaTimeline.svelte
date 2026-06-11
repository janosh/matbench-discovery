<script lang="ts">
  // Benchmark progress over time: each model plotted at its submission date, with a
  // step line tracing the running best ("SOTA frontier") of the selected metric.
  import { goto } from '$app/navigation'
  import type { ModelData } from '$lib'
  import { MODELS } from '$lib'
  import { get_nested_number, is_finite_num } from '$lib/metrics'
  import { sota_frontier_indices, sota_step_line } from '$lib/sota'
  import { ScatterPlot } from 'matterviz/plot'
  import type { DataSeries } from 'matterviz/plot'
  import type { HTMLAttributes } from 'svelte/elements'

  const METRIC_KEYS = [`CPS`, `F1`] as const
  type MetricKey = (typeof METRIC_KEYS)[number]
  interface TimelineMeta extends Record<string, unknown> {
    model_name: string
    model_key?: string
    date_added: string
    is_record: boolean
  }

  let {
    models = MODELS,
    model_filter = () => true,
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    models?: ModelData[]
    model_filter?: (model: ModelData) => boolean
  } = $props()

  let metric = $state<MetricKey>(`CPS`)

  const metric_configs: Record<
    MetricKey,
    { label: string; get: (model: ModelData) => number | null }
  > = {
    CPS: {
      label: `CPS (Combined Performance Score)`,
      get: (model) => model.CPS ?? null,
    },
    F1: {
      label: `F1 (discovery, unique prototypes)`,
      get: (model) =>
        get_nested_number(model, `metrics.discovery.unique_prototypes.F1`) ?? null,
    },
  }

  let points = $derived(
    models
      .filter(model_filter)
      .map((model) => ({
        model,
        date: Date.parse(model.date_added ?? ``),
        value: metric_configs[metric].get(model),
      }))
      .filter((pt): pt is typeof pt & { value: number } =>
        is_finite_num(pt.value) && is_finite_num(pt.date)
      ),
  )
  // sota_frontier_indices returns records in date order, so the step line can use
  // them directly without re-sorting
  let record_order = $derived(sota_frontier_indices(points, `higher`))
  let is_record = $derived.by(() => {
    const flags = points.map(() => false)
    for (const idx of record_order) flags[idx] = true
    return flags
  })
  let frontier = $derived(
    sota_step_line(record_order.map((idx) => points[idx]), Date.now()),
  )

  let series = $derived.by((): DataSeries<TimelineMeta>[] => {
    const model_points: DataSeries<TimelineMeta> = {
      x: points.map((pt) => pt.date),
      y: points.map((pt) => pt.value),
      label: `Models`,
      markers: `points`,
      metadata: points.map((pt, idx) => ({
        model_name: pt.model.model_name,
        model_key: pt.model.model_key,
        date_added: pt.model.date_added,
        is_record: is_record[idx],
      })),
      point_style: points.map((pt, idx) => ({
        fill: pt.model.color ?? `#4dabf7`,
        radius: is_record[idx] ? 6 : 4,
        stroke: is_record[idx] ? `white` : `transparent`,
        stroke_width: is_record[idx] ? 1 : 0,
      })),
      // label only frontier-defining models to keep the plot readable
      point_label: points.map((pt, idx) =>
        is_record[idx]
          ? { text: pt.model.model_name, font_size: `11px`, auto_placement: true }
          : {}
      ),
    }
    const frontier_line: DataSeries<TimelineMeta> = {
      x: frontier.x,
      y: frontier.y,
      label: `SOTA frontier`,
      markers: `line`,
      line_style: { stroke: `#4dabf7`, stroke_width: 2, line_dash: `6 3` },
    }
    return [frontier_line, model_points]
  })
</script>

<div class="metric-toggle">
  {#each METRIC_KEYS as key (key)}
    <button
      class:active={metric === key}
      onclick={() => (metric = key)}
      aria-pressed={metric === key}
    >
      {key}
    </button>
  {/each}
</div>

<ScatterPlot
  {series}
  x_axis={{ label: `Date added`, scale_type: `time`, format: `%b %Y` }}
  y_axis={{
    label: metric_configs[metric].label,
    format: `.2f`,
    range: [0, 1],
  }}
  legend={null}
  label_placement_config={{
    leader_line_threshold: 15,
    max_neighbors: { count: 3, radius: 40 },
  }}
  point_events={{
    onclick: ({ point }) => {
      const key = point.metadata?.model_key
      if (typeof key === `string`) goto(`/models/${encodeURIComponent(key)}`)
    },
  }}
  {...rest}
>
  {#snippet tooltip({ y_formatted, metadata })}
    {@const meta = metadata}
    {#if meta}
      <strong>{meta.model_name}</strong><br />
      Added: {meta.date_added}<br />
      {metric}: {y_formatted}
      {#if meta.is_record}<br /><em>set a new record</em>{/if}
    {:else}
      SOTA frontier: {y_formatted}
    {/if}
  {/snippet}
</ScatterPlot>

<style>
  .metric-toggle {
    display: flex;
    gap: 0.5em;
    justify-content: center;
    margin-bottom: 0.5em;
  }
  .metric-toggle button {
    padding: 2pt 8pt;
    font: inherit;
    border-radius: 4px;
    background: color-mix(in srgb, currentColor 8%, transparent);
  }
  .metric-toggle button.active {
    background: color-mix(in srgb, var(--link-color, #4dabf7) 25%, transparent);
  }
</style>
