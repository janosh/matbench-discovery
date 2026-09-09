<script module lang="ts">
  import type { SortState, UrlTableFilters } from '$lib/url-state.svelte'

  // the table's default sort; pages binding `sort` reuse this so URL sort params
  // are omitted when the table is at its resting state
  export const DEFAULT_TABLE_SORT: SortState = { column: `CPS`, dir: `desc` }

  // Shared HeatmapTable theme for striped sticky cells and flush-left row numbers.
  export const METRICS_TABLE_ROOT_STYLE = `--heatmap-sticky-cell-odd-bg: linear-gradient(var(--table-odd), var(--table-odd)), var(--page-bg); --heatmap-row-num-padding-left: 0;`
</script>

<script lang="ts">
  import { goto } from '$app/navigation'
  import OrgLogos from '$lib/model/OrgLogos.svelte'
  import TableControls from '$lib/table/TableControls.svelte'
  import {
    append_better_hint,
    metric_better_as,
    missing_metric_reason,
  } from '$lib/metrics'
  import {
    comparison,
    mark_compared_rows,
    row_model_key,
    toggle_row_model,
  } from '$lib/model-comparison.svelte'
  import { make_table_filters, MODELS } from '$lib/models.svelte'
  import type { DiscoverySet, Label, ModelData, SortDir } from '$lib/types'
  import type { CellSnippet, CellSnippetArgs, Column, RowData } from 'matterviz/table'
  import { HeatmapTable, is_invalid } from 'matterviz/table'
  import { format_num } from 'matterviz/labels'
  import { ActionMenu, type CmdAction, Icon } from 'svelte-widgets'
  import {
    Code,
    Download,
    Graph,
    Paper,
    PullRequest,
    Unavailable,
  } from 'svelte-widgets/icons'
  import { click_outside } from 'svelte-widgets/attachments'
  import type { HTMLAttributes } from 'svelte/elements'
  import { ALL_METRICS, HYPERPARAMS, METADATA_COLS } from '../labels'
  import { assemble_row_data } from '../metrics'

  type MetricsRow = ReturnType<typeof assemble_row_data>[number]
  type LinkData = MetricsRow[`Links`]
  type PredFilesDropdown = LinkData[`pred_files`] & { x: number; y: number }
  const resource_links = [
    [`paper`, `Read model paper`, Paper],
    [`repo`, `View source code`, Code],
    [`pr_url`, `View pull request`, PullRequest],
    [`checkpoint`, `Download model checkpoint`, Download],
  ] as const
  const format_element_count = (elements: string[]) =>
    `${elements.length}${elements.length ? ` (${elements.join(`, `)})` : ``}`
  const cells: Record<string, CellSnippet> = {
    ...Object.fromEntries(
      Object.values(ALL_METRICS).map(({ key }) => [key, metric_cell]),
    ),
    Links: links_cell,
    Org: affiliation_cell,
  }

  const { model_name, training_sets, targets, benchmark_added, links } = METADATA_COLS
  const { checkpoint_license, code_license, org } = METADATA_COLS
  const { graph_construction_radius, model_params } = HYPERPARAMS
  const heatmap_disabled_cols = new Set([
    training_sets.key,
    graph_construction_radius.key,
    benchmark_added.key,
    model_params.key,
  ])

  const default_columns = [
    model_name,
    ...Object.values(ALL_METRICS),
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

  let {
    discovery_set = $bindable(`unique_prototypes`),
    model_filter = $bindable(() => true),
    col_filter = $bindable(() => true),
    column_labels = default_columns,
    show_row_numbers = true,
    filters = make_table_filters(),
    column_order = $bindable([]),
    sort = $bindable({ ...DEFAULT_TABLE_SORT }),
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    discovery_set?: DiscoverySet
    model_filter?: (model: ModelData) => boolean
    col_filter?: (col: Label) => boolean
    column_labels?: Label[]
    show_row_numbers?: boolean
    filters?: UrlTableFilters
    column_order?: string[]
    sort?: { column: string; dir: SortDir }
  } = $props()
  // toggled from TableControls; no page binds it, so plain local state
  let show_selected_only = $state(false)

  let pred_files_dropdown = $state<PredFilesDropdown | null>(null)

  let metrics_data = $derived(
    mark_compared_rows(
      assemble_row_data(discovery_set, model_filter, filters.matches),
      show_selected_only,
    ),
  )
  let columns = $derived(
    column_labels.map((col): Column => {
      const better = col.better ?? metric_better_as(col.label) ?? undefined
      return {
        ...col,
        id: col.group ? `${col.key} (${col.group})` : col.key,
        cell: cells[col.key],
        ...(column_labels === default_columns && {
          color_scale: heatmap_disabled_cols.has(col.key) ? null : col.color_scale,
          ...(col === model_name && { style: `padding-left: 0;${col.style ?? ``}` }),
        }),
        better,
        description: append_better_hint(col, better),
        visible: col.visible !== false && col_filter(col),
      }
    }),
  )

  type ButtonMouseEvent = MouseEvent & { currentTarget: HTMLButtonElement }
  function show_dropdown(event: ButtonMouseEvent, link_data: LinkData) {
    event.stopPropagation()

    // position the dropdown at the button's document coordinates
    const rect = event.currentTarget.getBoundingClientRect()
    pred_files_dropdown = {
      ...link_data.pred_files,
      x: rect.left + globalThis.scrollX,
      y: rect.bottom + globalThis.scrollY,
    }
  }
  const close_dropdown = () => (pred_files_dropdown = null)

  let at = $state<{ x: number; y: number } | null>(null)
  let model_key = $state(``)
  let model = $derived(MODELS.find((md) => md.model_key === model_key))
  let selected = $derived(comparison.keys.has(model_key))

  let actions: CmdAction[] = $derived.by(() => {
    if (!model) return []
    const { model_name: name } = model
    const n_selected = comparison.keys.size
    return [
      {
        id: `toggle`,
        label: selected ? `Remove ${name} from comparison` : `Add ${name} to comparison`,
        action: () => comparison.toggle(model_key),
      },
      {
        id: `open`,
        label:
          selected && n_selected > 1
            ? `Compare ${n_selected} selected models`
            : `Compare ${name} with…`,
        action: () => comparison.open_with(model_key),
      },
      {
        id: `page`,
        label: `Open ${name} model page`,
        action: () => void goto(`/models/${model_key}`),
      },
    ]
  })

  function open_menu(event: MouseEvent) {
    // Links, buttons and inputs keep the browser's own menu (open in new tab, copy link).
    const target = event.target instanceof Element ? event.target : null
    if (!target || target.closest(`a, button, input, select`)) return
    const key = row_model_key(target.closest(`tbody tr`))
    if (!key) return
    event.preventDefault()
    event.stopPropagation() // capture phase: preempts HeatmapTable's column menu
    model_key = key
    at = { x: event.clientX, y: event.clientY }
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

{#snippet metric_cell({ row, col, val }: CellSnippetArgs)}
  {@const coverage =
    col.key === ALL_METRICS.pbe_vib_freq_error.key
      ? (row as MetricsRow).model.metrics?.diatomics?.pbe_vib_freq_coverage
      : undefined}
  {#if coverage}
    <span
      data-title={`Valid fits: ${coverage.n_valid}/${coverage.n_eligible} reference-eligible elements. No valid fit: ${format_element_count(coverage.failed_elements)}. Unavailable curves: ${format_element_count(coverage.missing_elements)}.`}
    >
      {typeof val === `number` && !is_invalid(val)
        ? format_num(val, col.format ?? `.3f`)
        : `n/a`}
      <small style="font-weight: 400; opacity: 0.75;"
        >· {coverage.n_valid}/{coverage.n_eligible}</small
      >
    </span>
  {:else if is_invalid(val)}
    <span data-title={missing_metric_reason((row as MetricsRow).model, col as Label)}
      >n/a</span
    >
  {:else if typeof val === `number`}
    {format_num(val, col.format ?? `.3f`)}
  {:else}
    {val}
  {/if}
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

<HeatmapTable
  data={metrics_data as RowData[]}
  row_key="model_key"
  {columns}
  bind:sort
  {show_row_numbers}
  default_num_format=".3f"
  bind:show_heatmap={filters.show_heatmap}
  bind:column_order
  export_data={{ formats: [`csv`], filename: `matbench-discovery-${discovery_set}` }}
  on_row_double_click={toggle_row_model}
  {...rest}
  oncontextmenucapture={open_menu}
  class={[`leaderboard`, rest.class]}
  root_style={METRICS_TABLE_ROOT_STYLE}
>
  {#snippet controls()}
    <TableControls bind:columns bind:show_selected_only {filters} />
  {/snippet}
</HeatmapTable>

<!-- Press dismissal prevents the right-click's mouseup from immediately closing the menu. -->
<ActionMenu
  {actions}
  bind:at
  trigger="none"
  dismiss={{ dismiss_on: `press` }}
  aria-label="Model row actions"
  style="font-size: 12px; --action-menu-padding: 0; --action-menu-item-padding: 3pt 8pt"
/>

{#if pred_files_dropdown}
  {@const { x, y, name, files } = pred_files_dropdown}
  {@const style = `position: absolute; left: ${x}px; top: ${y}px;`}
  <div
    class="pred-files-dropdown"
    {style}
    {@attach click_outside({ callback: close_dropdown })}
  >
    <h4 id="files-for">Files for {name}</h4>
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
