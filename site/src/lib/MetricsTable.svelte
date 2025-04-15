<script lang="ts">
  import { HeatmapTable, TableControls } from '$lib'
  import { DEFAULT_CPS_CONFIG } from '$lib/combined_perf_score'
  import { metric_better_as } from '$lib/metrics'
  import { pretty_num } from 'elementari'
  import { click_outside } from 'svelte-zoo/actions'
  import { ALL_METRICS, METADATA_COLS, calculate_metrics_data } from './metrics'
  import type {
    CombinedMetricConfig,
    DiscoverySet,
    HeatmapColumn,
    LinkData,
    ModelData,
    RowData,
  } from './types'

  interface Props {
    discovery_set?: DiscoverySet
    model_filter?: (model: ModelData) => boolean
    col_filter?: (col: HeatmapColumn) => boolean
    show_energy_only?: boolean
    show_noncompliant?: boolean
    config?: CombinedMetricConfig
    [key: string]: unknown
  }
  let {
    discovery_set = `unique_prototypes`,
    model_filter = () => true,
    col_filter = () => true,
    show_energy_only = false,
    show_noncompliant = false,
    config = DEFAULT_CPS_CONFIG,
    ...rest
  }: Props = $props()

  let active_files: { name: string; url: string }[] = $state([])
  let active_model_name = $state(``)
  let active_dropdown_pos = $state<{ x: number; y: number; name: string } | null>(null)
  const [compliant_clr, noncompliant_clr] = [`#4caf50`, `#4682b4`]

  let columns = $derived(
    [...ALL_METRICS, ...METADATA_COLS]
      .map((col) => {
        const better = col.better ?? metric_better_as(col.label)

        // append better=higher/lower to tooltip if applicable
        let tooltip = col.tooltip || ``
        if (better === `higher` || better === `lower`) {
          tooltip = tooltip ? `${tooltip} (${better}=better)` : `${better}=better`
        }
        const visible = col.visible !== false && col_filter(col)

        return { ...col, better, tooltip, visible } as HeatmapColumn
      })
      // Ensure Model column comes first
      .sort((col1, _col2) => (col1.label === `Model` ? -1 : 1)),
  )

  let metrics_data = $state<RowData[]>([])
  // recalculate metrics_data whenever filter settings, props, or metric_config change
  $effect(() => {
    metrics_data = calculate_metrics_data(
      discovery_set,
      model_filter,
      show_energy_only,
      show_noncompliant,
      config,
      compliant_clr,
      noncompliant_clr,
    )
  })

  function show_dropdown(event: MouseEvent, links: LinkData) {
    event.stopPropagation()

    // Get button position for dropdown placement
    const button = event.currentTarget as HTMLElement
    const rect = button.getBoundingClientRect()

    active_model_name = links.pred_files.name
    active_files = links.pred_files.files

    // Position dropdown relative to the viewport
    active_dropdown_pos = { x: rect.left, y: rect.bottom, name: links.pred_files.name }
  }
  const close_dropdown = () => (active_dropdown_pos = null)
</script>

<svelte:window
  onkeydown={(event) => {
    if (event.key === `Escape` && active_dropdown_pos) {
      close_dropdown()
      event.preventDefault()
    }
  }}
/>

<HeatmapTable
  data={metrics_data}
  {columns}
  initial_sort_column="CPS"
  initial_sort_direction="desc"
  sort_hint="Click on column headers to sort table rows"
  {...rest}
>
  {#snippet controls()}
    <div class="controls-container">
      <div class="controls-row">
        {#if show_noncompliant}
          {#each [[compliant_clr, `Compliant`], [noncompliant_clr, `Non-compliant`]] as [clr, label] (label)}
            <div class="legend-item">
              <span class="color-swatch" style="background-color: {clr};"></span>
              {label}
            </div>
          {/each}
        {/if}
        <TableControls bind:show_energy_only bind:show_noncompliant bind:columns />
      </div>
    </div>
  {/snippet}

  {#snippet cell({ col, val })}
    {#if col.label === `Links` && val && typeof val === `object` && `paper` in val}
      {@const links = val as LinkData}
      {#each Object.entries(links).filter(([key]) => key !== `pred_files`) as [key, link] (JSON.stringify(link))}
        {#if `url` in link && ![`missing`, `not available`, ``, null, undefined].includes(link.url)}
          <a href={link.url} target="_blank" rel="noopener noreferrer" title={link.title}>
            {@html link.icon}
          </a>
        {:else}
          <span title="{key} not available">
            <svg><use href="#icon-unavailable"></use></svg>
          </span>
        {/if}
      {/each}
      {#if links?.pred_files}
        <button
          class="pred-files-btn"
          aria-label="Download model prediction files"
          onclick={(event) => show_dropdown(event, links)}
        >
          <svg><use href="#icon-graph"></use></svg>
        </button>
      {/if}
    {:else if typeof val === `number` && col.format}
      {pretty_num(val, col.format)}
    {:else if val === undefined || val === null}
      n/a
    {:else}
      {@html val}
    {/if}
  {/snippet}
</HeatmapTable>

{#if active_dropdown_pos}
  {@const { x, y } = active_dropdown_pos}
  <div
    class="pred-files-dropdown"
    style="position: fixed; left: {x}px; top: {y}px;"
    use:click_outside={{ callback: close_dropdown }}
  >
    <h4>Files for {active_model_name}</h4>
    <ol>
      {#each active_files as { name, url } (url)}
        <li>
          <a href={url} target="_blank" rel="noopener noreferrer">
            {@html name}
          </a>
        </li>
      {/each}
    </ol>
  </div>
{/if}

<style>
  .pred-files-btn {
    background: none;
    padding: 0;
  }
  .pred-files-dropdown {
    transform: translateX(-100%);
    margin-left: 20px;
    background: var(--light-bg, white);
    color: var(--text-color, black);
    border-radius: 5px;
    padding: 0.75em;
  }
  .pred-files-dropdown h4 {
    margin: 0;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
  }
  .pred-files-dropdown ol {
    margin: 0;
    padding-left: 1em;
  }
  div.controls-row {
    display: flex;
    gap: 1em;
    flex-wrap: wrap;
    font-size: 0.85em;
  }
  div.legend-item {
    display: flex;
    place-items: center;
    gap: 0.3em;
  }
  span.color-swatch {
    width: 14px;
    height: 14px;
    border-radius: 2px;
  }
</style>
