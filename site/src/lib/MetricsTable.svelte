<script lang="ts">
  import { HeatmapTable, IconList, TableControls } from '$lib'
  import { metric_better_as } from '$lib/metrics'
  import type { Snippet } from 'svelte'
  import { click_outside } from 'svelte-zoo/actions'
  import { ALL_METRICS, HYPERPARAMS, METADATA_COLS } from './labels'
  import { assemble_row_data } from './metrics'
  import type { CellSnippetArgs, DiscoverySet, Label, LinkData, ModelData } from './types'

  interface Props {
    discovery_set?: DiscoverySet
    model_filter?: (model: ModelData) => boolean
    col_filter?: (col: Label) => boolean
    show_energy_only?: boolean
    show_non_compliant?: boolean
    show_heatmap?: boolean
    show_compliant?: boolean
    active_files?: { name: string; url: string }[]
    active_model_name?: string
    pred_files_dropdown_pos?: { x: number; y: number; name: string } | null
    [key: string]: unknown
  }
  let {
    discovery_set = $bindable(`unique_prototypes`),
    model_filter = $bindable(() => true),
    col_filter = $bindable(() => true),
    show_energy_only = $bindable(false),
    show_non_compliant = $bindable(true),
    show_heatmap = $bindable(true),
    show_compliant = $bindable(true),
    active_files = $bindable([]),
    active_model_name = $bindable(``),
    pred_files_dropdown_pos = $bindable(null),
    ...rest
  }: Props = $props()

  const { model_name, training_set, targets, date_added, links } = METADATA_COLS
  const { checkpoint_license, code_license, org } = METADATA_COLS
  const { graph_construction_radius, model_params } = HYPERPARAMS

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
      org,
    ]
      .map((col) => {
        const better = col.better ?? metric_better_as(col.label)

        // append better=higher/lower to tooltip if applicable
        let description = col.description ?? ``
        if (better === `higher` || better === `lower`) {
          description = description
            ? `${description} (${better}=better)`
            : `${better}=better`
        }
        const visible = col.visible !== false && col_filter(col)

        return { ...col, better, description, visible } as Label
      })
      // Ensure Model column comes first
      .sort((col1, _col2) => (col1.label === `Model` ? -1 : 1)),
  )

  let metrics_data = $derived(
    assemble_row_data(
      discovery_set,
      model_filter,
      show_energy_only,
      show_non_compliant,
      show_compliant,
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
    pred_files_dropdown_pos = {
      x: rect.left,
      y: rect.bottom,
      name: links.pred_files.name,
    }
  }
  const close_dropdown = () => (pred_files_dropdown_pos = null)
</script>

<svelte:window
  onkeydown={(event) => {
    if (event.key === `Escape` && pred_files_dropdown_pos) {
      close_dropdown()
      event.preventDefault()
    }
  }}
/>

{#snippet affiliation_cell({ row }: CellSnippetArgs)}
  <IconList icons={(row as ModelData).org_logos} />
{/snippet}

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
  special_cells={{
    Links: links_cell as unknown as Snippet<[CellSnippetArgs]>,
    Org: affiliation_cell as unknown as Snippet<[CellSnippetArgs]>,
  }}
  default_num_format=".3f"
  bind:show_heatmap
  {...rest}
>
  {#snippet controls()}
    <TableControls
      bind:show_energy_only
      bind:columns
      bind:show_heatmap
      bind:show_compliant
      bind:show_non_compliant
    />
  {/snippet}
</HeatmapTable>

{#if pred_files_dropdown_pos}
  {@const { x, y } = pred_files_dropdown_pos}
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
</style>
