<script lang="ts">
  import { goto } from '$app/navigation'
  import { format_num, Icon } from 'matterviz'
  import { ScatterPlot } from 'matterviz/plot'
  import type { ComponentProps } from 'svelte'
  import type { Label } from './types'

  type GitHubData = {
    name: string
    repo: string
    stars: number
    forks: number
    commits_last_year: number
    contributors: number
    model_key?: string // URL slug for model detail page
  }

  let {
    github_data = [],
    show_model_labels = true,
    color_scale = { type: `log` },
    ...rest
  }: {
    github_data?: GitHubData[]
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

  type GitHubMetricKey = keyof Omit<GitHubData, `name` | `repo` | `model_key`>

  let axes = $state({
    x: github_metrics.forks as Label,
    y: github_metrics.stars as Label,
    color_value: github_metrics.commits_last_year as Label,
    size_value: github_metrics.contributors as Label,
  })

  let is_fullscreen = $state(false)
  let container_el = $state<HTMLDivElement | null>(null)

  type PlotPoint = {
    x: number
    y: number
    color_value: number
    size_value: number
    metadata: { name: string; repo: string; model_key?: string }
  }
  let plot_data = $derived(
    github_data
      .map((item) => ({
        x: item[axes.x.key as GitHubMetricKey],
        y: item[axes.y.key as GitHubMetricKey],
        color_value: item[axes.color_value.key as GitHubMetricKey],
        size_value: item[axes.size_value.key as GitHubMetricKey],
        metadata: { name: item.name, repo: item.repo, model_key: item.model_key },
      }))
      .filter(({ x, y, color_value, size_value }) =>
        [x, y, color_value, size_value].every((val) => val != null)
      ) as PlotPoint[],
  )
  // O(1) lookup for tooltip by model name
  let plot_data_by_name = $derived(
    new Map(plot_data.map((d) => [d.metadata.name, d])),
  )

  const point_style = { fill: `#4dabf7` }
  let series = $derived({
    x: plot_data.map((d) => d.x) as number[],
    y: plot_data.map((d) => d.y) as number[],
    markers: `points` as const,
    point_style,
    metadata: plot_data.map((d) => d.metadata),
    color_values: plot_data.map((d) => Math.max(1, d.color_value)) as number[],
    size_values: plot_data.map((d) => d.size_value) as number[],
    point_label: (show_model_labels ? plot_data : []).map((d) => ({
      text: d.metadata.name,
      font_size: `14px`,
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
  <ScatterPlot
    style={`height: ${is_fullscreen ? `100%` : `600px`}`}
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
    color_bar={{ title: axes.color_value.label, tick_format: axes.color_value.format }}
    {color_scale}
    point_events={{
      onclick: ({ point }) => {
        const slug = point.metadata?.model_key ?? point.metadata?.name
        if (typeof slug === `string`) goto(`/models/${encodeURIComponent(slug)}`)
      },
    }}
    {...rest}
  >
    <button
      onclick={() => {
        if (is_fullscreen) document.exitFullscreen()
        else container_el?.requestFullscreen?.()
      }}
      aria-label="{is_fullscreen ? `Exit` : `Enter`} fullscreen"
      title="{is_fullscreen ? `Exit` : `Enter`} fullscreen"
    >
      <Icon icon="{is_fullscreen ? `Exit` : ``}Fullscreen" />
    </button>

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

<style>
  div.fullscreen {
    background: var(--page-bg);
    overflow: auto;
    padding: 2em;
    box-sizing: border-box;
  }
  button[title$='fullscreen'] {
    position: absolute;
    top: 1em;
    right: 3.3em;
    display: flex;
    padding: 8px;
    border-radius: 50%;
  }
</style>
