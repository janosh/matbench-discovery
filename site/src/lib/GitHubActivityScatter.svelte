<script lang="ts">
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
  }

  let { github_data = [], show_model_labels = true, ...rest }: {
    github_data?: GitHubData[]
    show_model_labels?: boolean
  } & ComponentProps<typeof ScatterPlot> = $props()

  const github_metrics = {
    stars: {
      key: `stars`,
      label: `GitHub Stars`,
      better: `higher`,
      format: `,.0f`,
      description: ``,
    },
    forks: {
      key: `forks`,
      label: `GitHub Forks`,
      better: `higher`,
      format: `,.0f`,
      description: ``,
    },
    commits_last_year: {
      key: `commits_last_year`,
      label: `Commits Last Year`,
      format: `,.0f`,
      description: ``,
    },
    contributors: {
      key: `contributors`,
      label: `Contributors`,
      better: `higher`,
      format: `,.0f`,
      description: ``,
    },
  } as const satisfies Record<string, Label>

  type GitHubMetricKey = keyof Omit<GitHubData, `name` | `repo`>

  let axes = $state({
    x: github_metrics.forks as Label,
    y: github_metrics.stars as Label,
    color_value: github_metrics.commits_last_year as Label,
    size_value: github_metrics.contributors as Label,
  })

  let is_fullscreen = $state(false)
  let container_el = $state<HTMLDivElement | null>(null)

  let plot_data = $derived(
    github_data
      .map((item) => ({
        x: item[axes.x.key as GitHubMetricKey],
        y: item[axes.y.key as GitHubMetricKey],
        color_value: item[axes.color_value.key as GitHubMetricKey],
        size_value: item[axes.size_value.key as GitHubMetricKey],
        metadata: { name: item.name, repo: item.repo },
      }))
      .filter(({ x, y, color_value, size_value }) =>
        [x, y, color_value, size_value].every((val) => val != null)
      ),
  )

  let series = $derived({
    x: plot_data.map((d) => d.x),
    y: plot_data.map((d) => d.y),
    markers: `points` as const,
    point_style: plot_data.map(() => ({ fill: `#4dabf7` })),
    metadata: plot_data.map((d) => d.metadata),
    color_values: plot_data.map((d) => d.color_value),
    size_values: plot_data.map((d) => d.size_value),
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
  style:height={is_fullscreen ? `100%` : `600px`}
  style="margin-block: 2em"
>
  <ScatterPlot
    series={[series]}
    x_axis={{ label: axes.x.label, format: axes.x.format, label_shift: { y: -50 } }}
    y_axis={{
      label: axes.y.label,
      format: axes.y.format,
      label_shift: { x: 50, y: -10 },
    }}
    color_bar={{ title: axes.color_value.label, tick_format: axes.color_value.format }}
    point_events={{
      onclick: ({ point }) =>
        point.metadata?.repo &&
        globalThis.open(`https://github.com/${point.metadata.repo}`, `_blank`),
    }}
    {...rest}
  >
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
        {@const point = plot_data.find((d) => d.metadata.name === metadata.name)}
        <strong>{metadata.name}</strong><br />
        {axes.x.label}: {x_formatted}<br />
        {axes.y.label}: {y_formatted}<br />
        {#if point}
          {axes.color_value.label}: {format_num(point.color_value)}<br />
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
  button {
    position: absolute;
    top: 20px;
    right: 3.5em;
    display: flex;
    padding: 8px;
    border-radius: 50%;
  }
</style>
