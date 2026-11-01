<script lang="ts">
  import { goto } from '$app/navigation'
  import { format_num } from 'matterviz/labels'
  import { ScatterPlot } from 'matterviz/plot'
  import type { ComponentProps } from 'svelte'
  import type { GitHubActivityData, Label } from '#lib/types.js'

  let {
    github_data,
    color_scale = { type: `log` },
    ...rest
  }: {
    github_data: GitHubActivityData[]
  } & ComponentProps<typeof ScatterPlot> = $props()

  const fmt = `,.0f`
  const github_metrics = {
    stars: { key: `stars`, label: `GitHub Stars`, format: fmt },
    forks: { key: `forks`, label: `GitHub Forks`, format: fmt },
    commits_last_year: {
      key: `commits_last_year`,
      label: `Commits Last Year`,
      format: fmt,
    },
    contributors: { key: `contributors`, label: `Contributors`, format: fmt },
  } as const satisfies Record<string, Partial<Label>>

  type GitHubMetric = (typeof github_metrics)[keyof typeof github_metrics]

  const axes: Record<`x` | `y` | `color_value` | `size_value`, GitHubMetric> = {
    x: github_metrics.forks,
    y: github_metrics.stars,
    color_value: github_metrics.commits_last_year,
    size_value: github_metrics.contributors,
  }

  let plot_data = $derived(
    github_data.flatMap((item) => {
      const point = {
        x: item[axes.x.key],
        y: item[axes.y.key],
        color_value: item[axes.color_value.key],
        size_value: item[axes.size_value.key],
        metadata: {
          name: item.name,
          repo: item.repo,
          model_key: item.model_key,
          commits_last_year: item.commits_last_year,
          contributors: item.contributors,
        },
      }
      const values = [point.x, point.y, point.color_value, point.size_value]
      if (!values.every((val) => typeof val === `number` && isFinite(val))) return []
      return [point]
    }),
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
    point_label: plot_data.map((point) => ({
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
        const slug = point.metadata?.model_key
        if (typeof slug === `string`) goto(`/models/${encodeURIComponent(slug)}`)
      },
    }}
    {...rest}
  >
    {#snippet tooltip({ x_formatted, y_formatted, metadata })}
      {#if metadata}
        {@const { name, commits_last_year, contributors } = metadata}
        <strong>{name}</strong><br />
        {axes.x.label}: {x_formatted}<br />
        {axes.y.label}: {y_formatted}<br />
        {#if typeof commits_last_year === `number`}
          {axes.color_value.label}: {format_num(commits_last_year)}<br />
        {/if}
        {#if typeof contributors === `number`}
          {axes.size_value.label}: {format_num(contributors)}<br />
        {/if}
      {/if}
    {/snippet}
  </ScatterPlot>
</div>
