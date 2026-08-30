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
  import { mark_compared_rows, toggle_row_model } from '$lib/model-comparison.svelte'
  import { make_table_filters } from '$lib/models.svelte'
  import type { UrlTableFilters } from '$lib/url-state.svelte'
  import type { DiscoverySet, Label, ModelData, SortDir, TableLabel } from '$lib/types'
  import type { CellSnippetArgs, Label as MattervizLabel, RowData } from 'matterviz'
  import { HeatmapTable } from 'matterviz'
  import { Icon } from 'svelte-widgets'
  import {
    Code,
    Download,
    Graph,
    Paper,
    PullRequest,
    Unavailable,
  } from 'svelte-widgets/icons'
  import { click_outside, tooltip } from 'svelte-widgets/attachments'
  import { untrack } from 'svelte'
  import type { HTMLAttributes } from 'svelte/elements'
  import { SvelteMap, SvelteSet } from 'svelte/reactivity'
  import { ALL_METRICS, HYPERPARAMS, METADATA_COLS } from '../labels'
  import { assemble_row_data } from '../metrics'
  import { heatmap_class } from '../table-export'

  type HeaderLabel = MattervizLabel & { tooltip_description?: string }
  type MetricsRow = ReturnType<typeof assemble_row_data>[number]
  type LinkData = MetricsRow[`Links`]
  type PredFilesDropdown = LinkData[`pred_files`] & { x: number; y: number }
  const resource_links = [
    [`paper`, `Read model paper`, Paper],
    [`repo`, `View source code`, Code],
    [`pr_url`, `View pull request`, PullRequest],
    [`checkpoint`, `Download model checkpoint`, Download],
  ] as const

  let {
    discovery_set = $bindable(`unique_prototypes`),
    model_filter = $bindable(() => true),
    col_filter = $bindable(() => true),
    filters = make_table_filters(),
    show_selected_only = $bindable(false),
    column_order = $bindable([]),
    sort = $bindable({ ...DEFAULT_TABLE_SORT }),
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    discovery_set?: DiscoverySet
    model_filter?: (model: ModelData) => boolean
    col_filter?: (col: Label) => boolean
    filters?: UrlTableFilters
    show_selected_only?: boolean
    column_order?: string[]
    sort?: { column: string; dir: SortDir }
  } = $props()

  const { model_name, training_sets, targets, benchmark_added, links } = METADATA_COLS
  const { checkpoint_license, code_license, org } = METADATA_COLS
  const { graph_construction_radius, model_params } = HYPERPARAMS
  const pinned_col_rank = (col: Label): number => (col.label === model_name.label ? 0 : 1)
  const heatmap_disabled_cols = new SvelteSet([
    training_sets.key,
    graph_construction_radius.key,
    benchmark_added.key,
    model_params.key,
  ])

  let pred_files_dropdown = $state<PredFilesDropdown | null>(null)

  // Reuse one row object per model across rebuilds: HeatmapTable keys its {#each}
  // by row-object identity, so its flip animation only runs when the SAME objects
  // reorder. Without this cache, every CPS/CMDS weight change rebuilds all rows and
  // re-sorts happen as delete+recreate with no row-movement animation. Cached rows
  // must be $state proxies: in-place updates on plain objects wouldn't trigger
  // fine-grained re-renders of changed cells (same object identity = no signal).
  const row_cache = new SvelteMap<string, MetricsRow>()
  function build_rows(): MetricsRow[] {
    const fresh_rows = mark_compared_rows(
      assemble_row_data(discovery_set, model_filter, filters.matches),
      show_selected_only,
    )
    // cache access is untracked so callers don't subscribe to the very row signals
    // this merge writes (which would re-trigger them and double-render the table)
    return untrack(() =>
      fresh_rows.map((row) => {
        const cached = row_cache.get(row.model_key)
        if (!cached) {
          const proxied = $state(row) // deep proxy for fine-grained cell updates
          row_cache.set(row.model_key, proxied)
          return proxied
        }
        // blank keys absent from the fresh row (e.g. after a discovery-set switch);
        // undefined renders/sorts like a missing key and avoids dynamic `delete`
        for (const key of Object.keys(cached)) {
          if (!(key in row)) Reflect.set(cached, key, undefined)
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
      benchmark_added,
      links,
      graph_construction_radius,
      checkpoint_license,
      code_license,
      training_sets,
      org,
    ]
      .map((col): TableLabel => {
        const better = col.better ?? metric_better_as(col.label) ?? undefined
        const visible = col.visible !== false && col_filter(col)
        return {
          ...col,
          better,
          color_scale: heatmap_disabled_cols.has(col.key) ? null : col.color_scale,
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
    columns.map((col): HeaderLabel => ({
      ...col,
      description: undefined,
      tooltip_description: col.description,
    })),
  )

  type ButtonMouseEvent = MouseEvent & { currentTarget: HTMLButtonElement }
  function show_dropdown(event: ButtonMouseEvent, link_data: LinkData) {
    event.stopPropagation()

    // Get button position for dropdown placement
    const rect = event.currentTarget.getBoundingClientRect()

    // Position dropdown relative to the button's position in the document
    pred_files_dropdown = {
      ...link_data.pred_files,
      x: rect.left + globalThis.scrollX,
      y: rect.bottom + globalThis.scrollY,
    }
  }
  const close_dropdown = () => (pred_files_dropdown = null)

  const header_tooltip = (content: string | undefined) => (node: Element) => {
    const header_cell = node.closest(`th`)
    return header_cell instanceof HTMLElement
      ? tooltip({ allow_html: true, content, placement: `top` })(header_cell)
      : undefined
  }
</script>

<svelte:window
  onkeydown={(event) => {
    if (event.key === `Escape` && pred_files_dropdown) {
      close_dropdown()
      event.preventDefault()
    }
  }}
/>

{#snippet affiliation_cell({ row }: CellSnippetArgs)}
  {@const metrics_row = row as MetricsRow}
  <OrgLogos org_logos={metrics_row.org_logos} authors={metrics_row.authors} />
{/snippet}

{#snippet links_cell({ val }: CellSnippetArgs)}
  {@const links = val as LinkData}
  {#each resource_links as [key, title, icon] (key)}
    {@const href = links[key]}
    {#if href}
      <a {href} target="_blank" rel="noopener noreferrer" {title}>
        <Icon {icon} />
      </a>
    {:else}
      <span title="{key} not available">
        <Icon icon={Unavailable} />
      </span>
    {/if}
  {/each}
  <button
    style="background: none; padding: 0"
    aria-label="Download model prediction files"
    onclick={(event) => show_dropdown(event, links)}
  >
    <Icon icon={Graph} />
  </button>
{/snippet}

{#snippet header_cell({ col }: { col: HeaderLabel })}
  <span
    class="header-label"
    style="display: inline-block"
    {@attach header_tooltip(col.tooltip_description)}
  >
    {@html col.label}
  </span>
{/snippet}

<HeatmapTable
  data={metrics_data as RowData[]}
  columns={table_columns}
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
  on_row_double_click={toggle_row_model}
  {...rest}
  root_style={METRICS_TABLE_ROOT_STYLE}
>
  {#snippet controls()}
    <TableControls bind:columns bind:show_selected_only {filters} />
  {/snippet}
</HeatmapTable>

{#if pred_files_dropdown}
  {@const { x, y, name, files } = pred_files_dropdown}
  {@const style = `position: absolute; left: ${x}px; top: ${y}px;`}
  <div
    class="pred-files-dropdown"
    {style}
    {@attach click_outside({ callback: close_dropdown })}
  >
    <h4>Files for {name}</h4>
    <ol>
      {#each files as { name: file_name, url } (url)}
        <li>
          <a href={url} target="_blank" rel="noopener noreferrer">
            {@html file_name}
          </a>
        </li>
      {/each}
    </ol>
  </div>
{/if}

<style>
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
