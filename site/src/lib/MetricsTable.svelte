<script lang="ts">
  import { HeatmapTable, TableControls } from '$lib'
  import { metric_better_as } from '$lib/metrics'
  import type { Snippet } from 'svelte'
  import { click_outside } from 'svelte-zoo/actions'
  import { ALL_METRICS, HYPERPARAMS, INFO_COLS, METADATA_COLS } from './labels'
  import { assemble_row_data } from './metrics'
  import type {
    CellSnippetArgs,
    DiscoverySet,
    LinkData,
    Metric,
    ModelData,
  } from './types'

  interface Props {
    discovery_set?: DiscoverySet
    model_filter?: (model: ModelData) => boolean
    col_filter?: (col: Metric) => boolean
    show_energy_only?: boolean
    show_noncompliant?: boolean
    [key: string]: unknown
  }
  let {
    discovery_set = `unique_prototypes`,
    model_filter = () => true,
    col_filter = () => true,
    show_energy_only = false,
    show_noncompliant = false,
    ...rest
  }: Props = $props()

  let active_files: { name: string; url: string }[] = $state([])
  let active_model_name = $state(``)
  let active_dropdown_pos = $state<{ x: number; y: number; name: string } | null>(null)
  const [compliant_clr, noncompliant_clr] = [`#4caf50`, `#4682b4`]
  const { model_name, training_set, model_params, targets, date_added, links } =
    METADATA_COLS
  const { graph_construction_radius } = HYPERPARAMS
  const { checkpoint_license, code_license } = INFO_COLS

  let columns = $derived(
    [
      ...Object.values(ALL_METRICS),
      model_name,
      training_set,
      model_params,
      targets,
      date_added,
      links,
      graph_construction_radius,
      checkpoint_license,
      code_license,
    ]
      .map((col) => {
        const better = col.better ?? metric_better_as(col.label)

        // append better=higher/lower to tooltip if applicable
        let description = col.description || ``
        if (better === `higher` || better === `lower`) {
          description = description
            ? `${description} (${better}=better)`
            : `${better}=better`
        }
        const visible = col.visible !== false && col_filter(col)

        return { ...col, better, description, visible } as Metric
      })
      // Ensure Model column comes first
      .sort((col1, _col2) => (col1.label === `Model` ? -1 : 1)),
  )

  let metrics_data = $derived(
    assemble_row_data(
      discovery_set,
      model_filter,
      show_energy_only,
      show_noncompliant,
      compliant_clr,
      noncompliant_clr,
    ),
  )

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

{#snippet links_cell({ val }: CellSnippetArgs)}
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
{/snippet}

<HeatmapTable
  data={metrics_data}
  {columns}
  initial_sort_column="CPS"
  initial_sort_direction="desc"
  sort_hint="Click on column headers to sort table rows"
  special_cells={{ Links: links_cell as unknown as Snippet<[CellSnippetArgs]> }}
  default_num_format=".3f"
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
