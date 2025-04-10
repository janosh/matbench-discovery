<script lang="ts">
  import { HeatmapTable, TableControls, get_metric_rank_order } from '$lib'
  import { DEFAULT_CPS_CONFIG } from '$lib/combined_perf_score'
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
  let pred_file_modal: HTMLDialogElement | null = $state(null)
  const [compliant_clr, noncompliant_clr] = [`#4caf50`, `#4682b4`]

  let columns = $derived(
    [...ALL_METRICS, ...METADATA_COLS]
      .map((col) => {
        const better = col.better ?? get_metric_rank_order(col.label)

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
</script>

<svelte:window
  onkeydown={(event) => {
    if (event.key === `Escape` && pred_file_modal?.open) {
      pred_file_modal.open = false
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
      {#each Object.entries(links).filter(([key]) => key !== `pred_files`) as [key, link] (key)}
        {#if `url` in link && ![`missing`, `not available`, ``, null, undefined].includes(link.url)}
          <a href={link.url} target="_blank" rel="noopener noreferrer" title={link.title}>
            {link.icon}
          </a>
        {:else}
          <span title="{key} not available">🚫</span>
        {/if}
      {/each}
      {#if links?.pred_files}
        <button
          class="pred-files-btn"
          title="Download model prediction files"
          onclick={() => {
            if (!pred_file_modal) return
            pred_file_modal.open = true
            active_files = links.pred_files.files
            active_model_name = links.pred_files.name
          }}
        >
          📊
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

<dialog
  bind:this={pred_file_modal}
  use:click_outside={{
    callback: () => {
      if (pred_file_modal?.open) pred_file_modal.open = false
    },
  }}
>
  <div class="modal-content">
    <button
      class="close-btn"
      onclick={() => {
        if (pred_file_modal?.open) pred_file_modal.open = false
      }}
      title="Close (or click escape)"
    >
      ×
    </button>
    <h3>Download prediction files for {active_model_name}</h3>
    <ol class="pred-files-list">
      {#each active_files as { name, url } (name + url)}
        <li>
          <a href={url} target="_blank" rel="noopener noreferrer">
            {name}
          </a>
        </li>
      {/each}
    </ol>
  </div>
</dialog>

<style>
  dialog {
    visibility: hidden;
    opacity: 0;
    background: var(--light-bg);
    color: var(--text-color);
    border: none;
    border-radius: 5pt;
    padding: 0;
    max-width: min(90vw, 500px);
  }

  dialog[open] {
    visibility: visible;
    opacity: 1;
    z-index: 2;
  }

  .pred-files-btn {
    background: none;
    padding: 0;
  }

  .modal-content {
    padding: 1em;
  }

  .modal-content h3 {
    margin: 0 0 1ex;
  }

  .pred-files-list {
    margin: 0;
    padding: 0 1em;
  }

  .close-btn {
    position: absolute;
    top: 0;
    right: 0;
    background: none;
    cursor: pointer;
    font-size: 24px;
  }
  .close-btn:hover {
    color: var(--link-color);
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
