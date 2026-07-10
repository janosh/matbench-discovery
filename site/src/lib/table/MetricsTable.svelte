<script module lang="ts">
  import type { SortState } from '$lib/url-state.svelte'

  // the table's default sort; pages binding `sort` reuse this so URL sort params
  // are omitted when the table is at its resting state
  export const DEFAULT_TABLE_SORT: SortState = { column: `CPS`, dir: `desc` }

  // Shared HeatmapTable theme for striped sticky cells and flush-left row numbers.
  export const METRICS_TABLE_ROOT_STYLE = `--heatmap-sticky-cell-odd-bg: linear-gradient(var(--table-odd), var(--table-odd)), var(--page-bg); --heatmap-row-num-padding-left: 0;`
</script>

<script lang="ts">
  import { OrgLogos, TableControls } from '$lib'
  import { append_better_hint, metric_better_as } from '$lib/metrics'
  import { make_table_filters } from '$lib/models.svelte'
  import type { UrlTableFilters } from '$lib/url-state.svelte'
  import type { DiscoverySet, Label, LinkData, ModelData, SortDir } from '$lib/types'
  import type { CellSnippetArgs, Label as MattervizLabel } from 'matterviz'
  import { HeatmapTable, Icon } from 'matterviz'
  import { click_outside, tooltip } from 'svelte-multiselect/attachments'
  import { untrack } from 'svelte'
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
    filters = make_table_filters(),
    show_selected_only = $bindable(false),
    active_files = $bindable([]),
    active_model_name = $bindable(``),
    pred_files_dropdown_pos = $bindable(null),
    selected_models = $bindable(new SvelteSet<string>()),
    column_order = $bindable([]),
    sort = $bindable({ ...DEFAULT_TABLE_SORT }),
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    discovery_set?: DiscoverySet
    model_filter?: (model: ModelData) => boolean
    col_filter?: (col: Label) => boolean
    filters?: UrlTableFilters
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

  // Reuse one row object per model across rebuilds: HeatmapTable keys its {#each}
  // by row-object identity, so its flip animation only runs when the SAME objects
  // reorder. Without this cache, every CPS/CMDS weight change rebuilds all rows and
  // re-sorts happen as delete+recreate with no row-movement animation. Cached rows
  // must be $state proxies: in-place updates on plain objects wouldn't trigger
  // fine-grained re-renders of changed cells (same object identity = no signal).
  type MetricsRow = ReturnType<typeof assemble_row_data>[number]
  const row_cache = new Map<string, MetricsRow>()
  function build_rows(): MetricsRow[] {
    // tracked snapshot of selections: reads inside the untrack block below wouldn't
    // subscribe, but selection changes must re-run this sync
    const selected_names = new Set(selected_models)
    const fresh_rows = assemble_row_data(
      discovery_set,
      model_filter,
      filters.matches,
    ).filter((row) => !show_selected_only || selected_names.has(row.model_name))
    // cache access is untracked so callers don't subscribe to the very row signals
    // this merge writes (which would re-trigger them and double-render the table)
    return untrack(() =>
      fresh_rows.map((row) => {
        // Only apply selected styles when not filtering to show only selected models
        row.class =
          !show_selected_only && selected_names.has(row.model_name) ? `highlight` : null
        const cached = row_cache.get(row.model_name)
        if (!cached) {
          const proxied = $state(row) // deep proxy for fine-grained cell updates
          row_cache.set(row.model_name, proxied)
          return proxied
        }
        // blank keys absent from the fresh row (e.g. after a discovery-set switch);
        // undefined renders/sorts like a missing key and avoids dynamic `delete`
        for (const key of Object.keys(cached)) {
          if (!(key in row)) (cached as Record<string, unknown>)[key] = undefined
        }
        return Object.assign(cached, row)
      }),
    )
  }
  // initialized eagerly (so SSR/prerendered HTML isn't empty), then synced by
  // $effect.pre. NOT a $derived: the cache merge writes $state proxies, and state
  // writes during derived evaluation schedule a second render flush -- the table then
  // reconciled twice per weight change and the second pass cancelled the first's flip
  // animations (rows jumped instantly). $effect.pre settles everything in one flush.
  let metrics_data = $state(build_rows())
  $effect.pre(() => {
    metrics_data = build_rows()
  })
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
          // tuck the Model cells (always adjacent to the rank column, being pinned
          // first) against the rank numbers
          ...(col === model_name && { style: `padding-left: 0;${col.style ?? ``}` }),
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
  {@const { org_logos = [], authors = [] } = row as Pick<
    ModelData,
    `org_logos` | `authors`
  >}
  <OrgLogos {org_logos} {authors} />
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
  show_row_numbers
  default_num_format=".3f"
  bind:show_heatmap={filters.show_heatmap}
  bind:column_order
  {heatmap_class}
  {header_cell}
  onrowdblclick={(event, row) => {
    if (typeof row.model_name === `string`) {
      handle_row_double_click(event, row.model_name)
    }
  }}
  {...rest}
  root_style={METRICS_TABLE_ROOT_STYLE}
>
  {#snippet controls()}
    <TableControls bind:columns bind:show_selected_only {filters} {selected_count} />
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
