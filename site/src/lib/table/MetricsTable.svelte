<script lang="ts">
  import { OrgLogos, TableControls } from '$lib'
  import { append_better_hint, metric_better_as } from '$lib/metrics'
  import type { DiscoverySet, Label, LinkData, ModelData, SortDir } from '$lib/types'
  import type { CellSnippetArgs, Label as MattervizLabel } from 'matterviz'
  import { HeatmapTable, Icon } from 'matterviz'
  import { click_outside, tooltip } from 'svelte-multiselect/attachments'
  import type { HTMLAttributes } from 'svelte/elements'
  import { SvelteSet } from 'svelte/reactivity'
  import { ALL_METRICS, HYPERPARAMS, METADATA_COLS } from '../labels'
  import { assemble_row_data } from '../metrics'
  import { heatmap_class } from '../table-export'

  type HeaderLabel = MattervizLabel & { tooltip_description?: string }

  let {
    discovery_set = $bindable(`unique_prototypes`),
    model_filter = $bindable(() => true),
    col_filter = $bindable(() => true),
    show_energy_only = $bindable(false),
    show_non_compliant = $bindable(true),
    show_heatmap = $bindable(true),
    show_compliant = $bindable(true),
    show_selected_only = $bindable(false),
    active_files = $bindable([]),
    active_model_name = $bindable(``),
    pred_files_dropdown_pos = $bindable(null),
    selected_models = $bindable(new SvelteSet<string>()),
    column_order = $bindable([]),
    sort = $bindable({ column: `CPS`, dir: `desc` }),
    show_energy_only_toggle = false,
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    discovery_set?: DiscoverySet
    model_filter?: (model: ModelData) => boolean
    col_filter?: (col: Label) => boolean
    show_energy_only?: boolean
    show_energy_only_toggle?: boolean
    show_non_compliant?: boolean
    show_heatmap?: boolean
    show_compliant?: boolean
    show_selected_only?: boolean
    active_files?: { name: string; url: string }[]
    active_model_name?: string
    pred_files_dropdown_pos?: { x: number; y: number; name: string } | null
    selected_models?: Set<string>
    column_order?: string[]
    sort?: { column: string; dir: SortDir }
  } = $props()

  const { model_name, training_set, targets, date_added, links } = METADATA_COLS
  const { checkpoint_license, code_license, org } = METADATA_COLS
  const { graph_construction_radius, model_params } = HYPERPARAMS
  const pinned_col_rank = (col: Label): number => (col.label === model_name.label ? 0 : 1)

  let selected_count = $derived(selected_models.size)

  let metrics_data = $derived(
    assemble_row_data(
      discovery_set,
      model_filter,
      show_energy_only,
      show_non_compliant,
      show_compliant,
    )
      .filter((row) => {
        // If show_selected_only is true, only show selected models
        if (show_selected_only) {
          return selected_models.has(row.model_name)
        }
        return true
      })
      .map((row) => {
        const is_selected = selected_models.has(row.model_name)
        // Only apply selected styles when not filtering to show only selected models
        row.class = is_selected && !show_selected_only ? `highlight` : null
        return row
      }),
  )
  let columns = $derived(
    [
      ...Object.values(ALL_METRICS),
      model_name,
      model_params,
      targets,
      date_added,
      links,
      graph_construction_radius,
      checkpoint_license,
      code_license,
      training_set,
      org,
    ]
      .map((col): Label => {
        const better = col.better ?? metric_better_as(col.label) ?? undefined
        const visible = col.visible !== false && col_filter(col)
        return {
          ...col,
          better,
          description: append_better_hint(col, better),
          visible,
        }
      })
      // Keep the sticky model column first, preserving definition order for the rest.
      .toSorted((col1, col2) => pinned_col_rank(col1) - pinned_col_rank(col2)),
  )
  let table_columns = $derived(
    columns.map(
      (col): HeaderLabel => ({
        ...col,
        better: col.better ?? undefined,
        description: undefined,
        tooltip_description: col.description,
      }),
    ),
  )

  function show_dropdown(event: MouseEvent, link_data: LinkData) {
    event.stopPropagation()

    // Get button position for dropdown placement
    const button = event.currentTarget as HTMLElement
    const rect = button.getBoundingClientRect()

    active_model_name = link_data.pred_files.name
    active_files = link_data.pred_files.files

    // Position dropdown relative to the button's position in the document
    const { name } = link_data.pred_files
    pred_files_dropdown_pos = {
      x: rect.left + window.scrollX,
      y: rect.bottom + window.scrollY,
      name,
    }
  }
  const close_dropdown = () => (pred_files_dropdown_pos = null)

  function toggle_model_selection(row_model_name: string) {
    const new_selected = new SvelteSet(selected_models)
    if (new_selected.has(row_model_name)) new_selected.delete(row_model_name)
    else new_selected.add(row_model_name)
    selected_models = new_selected
  }

  function handle_row_double_click(event: MouseEvent, row_model_name: string) {
    event.preventDefault()
    toggle_model_selection(row_model_name)
  }

  const header_tooltip = (content: string | undefined) => (node: Element) => {
    const header_cell = node.closest(`th`)
    return header_cell instanceof HTMLElement
      ? tooltip({ allow_html: true, content, placement: `top` })(header_cell)
      : undefined
  }
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
  {@const model = row as ModelData}
  <OrgLogos org_logos={model.org_logos ?? []} authors={model.authors ?? []} />
{/snippet}

{#snippet links_cell({ val }: CellSnippetArgs)}
  {@const links = val as unknown as LinkData}
  {#if links}
    {#each Object.entries(links).filter(([key]) => key !== `pred_files`) as [key, link] (JSON.stringify(link))}
      {#if `url` in link && ![`missing`, `not available`, ``, null, undefined].includes(link.url)}
        <a href={link.url} target="_blank" rel="noopener noreferrer" title={link.title}>
          <Icon icon={link.icon} />
        </a>
      {:else}
        <span title="{key} not available">
          <Icon icon="Unavailable" />
        </span>
      {/if}
    {/each}
    {#if links.pred_files}
      <button
        class="pred-files-btn"
        aria-label="Download model prediction files"
        onclick={(event) => show_dropdown(event, links)}
      >
        <Icon icon="Graph" />
      </button>
    {/if}
  {/if}
{/snippet}

{#snippet header_cell({ col }: { col: HeaderLabel })}
  <span class="header-label" {@attach header_tooltip(col.tooltip_description)}>
    {@html col.label}
  </span>
{/snippet}

<HeatmapTable
  data={metrics_data}
  columns={table_columns as MattervizLabel[]}
  bind:sort
  special_cells={{
    Links: links_cell,
    Org: affiliation_cell,
  }}
  default_num_format=".3f"
  bind:show_heatmap
  bind:column_order
  {heatmap_class}
  {header_cell}
  onrowdblclick={(event, row) => {
    if (typeof row.model_name === `string`) {
      handle_row_double_click(event, row.model_name)
    }
  }}
  {...rest}
>
  {#snippet controls()}
    <TableControls
      bind:columns
      bind:show_energy_only
      bind:show_heatmap
      bind:show_compliant
      bind:show_non_compliant
      bind:show_selected_only
      {show_energy_only_toggle}
      {selected_count}
    />
  {/snippet}
</HeatmapTable>

{#if pred_files_dropdown_pos}
  {@const { x, y } = pred_files_dropdown_pos}
  {@const style = `position: absolute; left: ${x}px; top: ${y}px;`}
  <div
    class="pred-files-dropdown"
    {style}
    {@attach click_outside({ callback: close_dropdown })}
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
  .header-label {
    display: inline-block;
  }
  .pred-files-btn {
    background: none;
    padding: 0;
  }
  .pred-files-dropdown {
    transform: translateX(-100%);
    margin-left: 20px;
    background: var(--page-bg);
    border: 1px solid var(--border);
    border-radius: 5px;
    padding: 4pt 11pt;
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
