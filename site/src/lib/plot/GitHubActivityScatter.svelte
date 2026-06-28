<script lang="ts">
  import { goto } from '$app/navigation'
  import { format_num } from 'matterviz'
  import { ScatterPlot } from 'matterviz/plot'
  import type { ComponentProps } from 'svelte'
  import type { GitHubActivityData, Label } from '$lib/types'

  let {
    github_data = [],
    show_model_labels = true,
    color_scale = { type: `log` },
    ...rest
  }: {
    github_data?: GitHubActivityData[]
    show_model_labels?: boolean
  } & ComponentProps<typeof ScatterPlot> = $props()

  const fmt = `,.0f`
  const github_metrics = {
    stars: { key: `stars`, label: `GitHub Stars`, better: `higher`, format: fmt },
    forks: { key: `forks`, label: `GitHub Forks`, better: `higher`, format: fmt },
    commits_last_year: {
      key: `commits_last_year`,
      label: `Commits Last Year`,
      format: fmt,
    },
    contributors: {
      key: `contributors`,
      label: `Contributors`,
      better: `higher`,
      format: fmt,
    },
  } as const satisfies Record<string, Partial<Label>>

  type GitHubMetric = (typeof github_metrics)[keyof typeof github_metrics]

  let axes = $state<Record<`x` | `y` | `color_value` | `size_value`, GitHubMetric>>({
    x: github_metrics.forks,
    y: github_metrics.stars,
    color_value: github_metrics.commits_last_year,
    size_value: github_metrics.contributors,
  })

  let plot_data = $derived(
    github_data.flatMap((item) => {
      const point = {
        x: item[axes.x.key],
        y: item[axes.y.key],
        color_value: item[axes.color_value.key],
        size_value: item[axes.size_value.key],
        metadata: { name: item.name, repo: item.repo, model_key: item.model_key },
      }
      const values = [point.x, point.y, point.color_value, point.size_value]
      if (!values.every((val) => typeof val === `number` && isFinite(val))) return []
      return [point]
    }),
  )
  // O(1) lookup for tooltip by model name
  let plot_data_by_name = $derived(
    new Map(plot_data.map((point) => [point.metadata.name, point])),
  )

  const point_style = { fill: `#4dabf7` }
  let series = $derived({
    x: plot_data.map((point) => point.x),
    y: plot_data.map((point) => point.y),
    markers: `points` as const,
    point_style,
    metadata: plot_data.map((point) => point.metadata),
    color_values: plot_data.map((point) => Math.max(1, point.color_value)),
    size_values: plot_data.map((point) => point.size_value),
    point_label: (show_model_labels ? plot_data : []).map((point) => ({
      text: point.metadata.name,
      font_size: `11px`,
      auto_placement: true,
    })),
  })
</script>

<div class="bleed-1400" style="margin-block: 2em">
  <ScatterPlot
    style="height: 600px"
    series={[series]}
    x_axis={{
      label: axes.x.label,
      format: axes.x.format,
      label_shift: { y: -50 },
      range: [0, null],
    }}
    y_axis={{
      label: axes.y.label,
      format: axes.y.format,
      range: [0, null],
    }}
    color_bar={{ title: axes.color_value.label, tick_format: `~s` }}
    {color_scale}
    label_placement_config={{ max_neighbors: { count: 3, radius: 40 } }}
    point_events={{
      onclick: ({ point }) => {
        const slug = point.metadata?.model_key ?? point.metadata?.name
        if (typeof slug === `string`) goto(`/models/${encodeURIComponent(slug)}`)
      },
    }}
    {...rest}
  >
    {#snippet tooltip({ x_formatted, y_formatted, metadata })}
      {@const name = (metadata as { name?: string })?.name}
      {#if name}
        {@const point = plot_data_by_name.get(name)}
        <strong>{name}</strong><br />
        {axes.x.label}: {x_formatted}<br />
        {axes.y.label}: {y_formatted}<br />
        {#if point?.color_value != null}
          {axes.color_value.label}: {format_num(point.color_value)}<br />
        {/if}
        {#if point?.size_value != null}
          {axes.size_value.label}: {format_num(point.size_value)}<br />
        {/if}
      {/if}
    {/snippet}
  </ScatterPlot>
</div>
