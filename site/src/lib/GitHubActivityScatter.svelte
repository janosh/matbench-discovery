<script lang="ts">
  import { extent } from 'd3-array'
  import { DraggablePane, format_num, Icon } from 'matterviz'
  import type { D3ColorSchemeName } from 'matterviz/colors'
  import { ColorScaleSelect, ScatterPlot } from 'matterviz/plot'
  import Select from 'svelte-multiselect'
  import type { HTMLAttributes } from 'svelte/elements'
  import type { Label } from './types'

  type GitHubData = {
    name: string
    repo: string
    stars: number
    forks: number
    commits_last_year: number
    contributors: number
  }

  let {
    github_data = [],
    show_model_labels = true,
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    github_data?: GitHubData[]
    show_model_labels?: boolean
  } = $props()

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

  const options = Object.values(github_metrics)

  let axes = $state({
    x: github_metrics.forks,
    y: github_metrics.stars,
    color_value: github_metrics.commits_last_year,
    size_value: github_metrics.contributors,
  })
  let log = $state({ x: false, y: false, color_value: false, size_value: false })
  let is_fullscreen = $state(false)
  let show_extra_controls = $state(false)
  let container_el = $state<HTMLDivElement | null>(null)
  let color_scheme = $state<D3ColorSchemeName>(`Viridis`)
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

  let plot_data = $derived(
    github_data
      .map((item) => ({
        x: item[axes.x.key as keyof GitHubData] as number,
        y: item[axes.y.key as keyof GitHubData] as number,
        color_value: item[axes.color_value.key as keyof GitHubData] as number,
        size_value: item[axes.size_value.key as keyof GitHubData] as number,
        metadata: { name: item.name, repo: item.repo },
      }))
      .filter(
        ({ x, y, color_value, size_value }) =>
          x != null && y != null && color_value != null && size_value != null,
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
      offset_y: 0,
      offset_x: 10,
      font_size: `${label_font_size}px`,
      color: `black`,
      auto_placement: true,
    })),
  })
</script>

{#if github_data.length === 0}
  <div class="no-data">
    <p>GitHub activity data is currently unavailable.</p>
    <p style="font-size: 0.9em; color: gray">
      Data will be fetched during the next site build.
    </p>
  </div>
{:else}
  <div class="plot-container" bind:this={container_el}>
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

    <DraggablePane
      bind:show={show_extra_controls}
      toggle_props={{
        class: `settings-toggle`,
        style:
          `position: absolute; top: 20px; right: 20px; background-color: var(--btn-bg); border-radius: 50%; padding: 4pt;`,
      }}
      pane_props={{
        style: `border: 1px solid var(--border);`,
        'aria-hidden': !show_extra_controls,
      }}
    >
      <div style="display: grid; grid-template-columns: auto 1fr; gap: 8pt 1em">
        <label
          style="grid-column: 1/-1"
          title="Toggle visibility of repository name labels on the scatter plot points"
        >
          <input type="checkbox" bind:checked={show_model_labels} /> Show Labels
        </label>
        <ColorScaleSelect
          bind:value={color_scheme}
          selected={color_scheme}
          style="margin: 0; grid-column: 1/-1"
        />
        <label title="Toggle the visibility of vertical grid lines">
          <input type="checkbox" bind:checked={x_grid} /> X Grid
        </label>
        <label title="Set the approximate number of ticks on the X axis">
          Ticks:
          <input
            type="number"
            min="0"
            max="20"
            bind:value={x_ticks}
            style="width: 50px"
          />
        </label>

        <label title="Toggle the visibility of horizontal grid lines">
          <input type="checkbox" bind:checked={y_grid} /> Y Grid
        </label>
        <label title="Set the approximate number of ticks on the Y axis">
          Ticks:
          <input
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
        >
          Point Size
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
          title="Adjust the font size of the repository name labels (in pixels)"
        >
          Label Size
        </label>
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
        >
          Label Link
        </label>
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

    <div class="controls-grid">
      {#each [`x`, `y`, `color_value`, `size_value`] as const as control (control)}
        {@const control_label = control === `color_value`
        ? `Color`
        : control === `size_value`
        ? `Size`
        : control === `x`
        ? `X Axis`
        : `Y Axis`}
        {@const [min, max] = extent(plot_data, (d) => d[control])}
        {@const disable_log = min == null || max == null || min <= 0 || 100 * min > max}
        <label for={control}>{control_label}</label>
        <Select
          {options}
          id={control}
          selected={[axes[control]]}
          bind:value={axes[control]}
          placeholder="Select {control_label}"
          maxSelect={1}
          minSelect={1}
          style="width: 100%; max-width: none; margin: 0"
          liSelectedStyle="font-size: 16px;"
        >
          {#snippet children({ option })}
            {
              typeof option === `object` && `label` in option
              ? option.label
              : String(option)
            }
          {/snippet}
        </Select>
        <label
          aria-disabled={disable_log}
          style:visibility={disable_log ? `hidden` : `visible`}
        >
          <input type="checkbox" bind:checked={log[control]} disabled={disable_log} />
          Log scale
        </label>
      {/each}
    </div>

    <div
      class="{is_fullscreen ? `` : `bleed-1400`} {rest.class ?? ``}"
      style:height={is_fullscreen ? `100%` : `600px`}
      style="margin-block: 2em"
    >
      <ScatterPlot
        series={[series]}
        x_axis={{
          label: axes.x.label,
          format: axes.x.format,
          scale_type: log.x ? `log` : `linear`,
          label_shift: { y: -50 },
          ticks: x_ticks,
        }}
        y_axis={{
          label: axes.y.label,
          format: axes.y.format,
          scale_type: log.y ? `log` : `linear`,
          label_shift: { x: 50, y: -10 },
          ticks: y_ticks,
        }}
        display={{ x_grid, y_grid }}
        color_scale={{ scheme: color_scheme, type: log.color_value ? `log` : `linear` }}
        size_scale={{
          radius_range: [5 * size_multiplier, 20 * size_multiplier],
          type: log.size_value ? `log` : `linear`,
        }}
        color_bar={{
          title: `${axes.color_value.label}${
            axes.color_value.better ? ` (${axes.color_value.better}=better)` : ``
          }`,
          margin: { t: 30, l: 80, b: 80, r: 50 },
          tick_format: axes.color_value.format,
        }}
        label_placement_config={{
          link_strength,
          link_distance,
          link_distance_range: [min_link_distance, max_link_distance],
        }}
        point_events={{
          onclick: ({ point }) =>
            globalThis.open(
              `https://github.com/${point.metadata?.repo ?? ``}`,
              `_blank`,
            ),
        }}
        controls={{
          toggle_props: { style: `position: absolute; top: 10px; right: 50px` },
        }}
        {...rest}
      >
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
  </div>
{/if}

<style>
  div.no-data {
    padding: 2em;
    text-align: center;
    background-color: var(--card-bg);
    border-radius: 4px;
    margin-block: 2em;
  }
  div.plot-container {
    position: relative;
    --scatter-min-height: 600px;
  }
  div.plot-container:fullscreen {
    background: var(--page-bg);
    overflow: auto;
    padding: 2em;
    box-sizing: border-box;
    display: flex;
    flex-direction: column;
    gap: 1em;
  }
  button {
    position: absolute;
    top: 20px;
    right: -20px;
    display: flex;
    padding: 8px;
    border-radius: 50%;
  }
  button[title='Exit fullscreen'] {
    right: 3.5em;
  }
  div.controls-grid {
    display: grid;
    grid-template-columns: auto 1fr auto;
    gap: 1ex;
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
