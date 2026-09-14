<script module lang="ts">
  import type { SortState, UrlTableFilters } from '$lib/url-state.svelte'

  // Omit the default sort from shared URLs.
  export const DEFAULT_TABLE_SORT: SortState = { column: `CPS`, dir: `desc` }

  // Shared HeatmapTable theme for striped sticky cells and flush-left row numbers.
  export const METRICS_TABLE_ROOT_STYLE = `--heatmap-sticky-cell-odd-bg: linear-gradient(var(--table-odd), var(--table-odd)), var(--page-bg); --heatmap-row-num-padding-left: 0; --heatmap-column-max-width: 14.4em;`
</script>

<script lang="ts">
  import { goto } from '$app/navigation'
  import data_files from '$pkg/data-files.yml'
  import {
    bind_url_params,
    sort_from_query,
    sort_url_entries,
  } from '$lib/url-state.svelte'
  import { MediaQuery } from 'svelte/reactivity'
  import { onMount, untrack } from 'svelte'
  import { CPS_CONFIG } from '$lib/combined-scores.svelte'
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
  import {
    ACTIVE_MODELS,
    make_table_filters,
    MODELS,
    score_weight_records,
  } from '$lib/models.svelte'
  import type { DiscoverySet, Label, ModelData } from '$lib/types'
  import type { CellSnippet, CellSnippetArgs, Column, RowData } from 'matterviz/table'
  import {
    cell_text,
    HeatmapTable,
    is_invalid,
    sort_table_rows,
    table_to_delimited,
  } from 'matterviz/table'
  import { strip_html } from 'matterviz/utils'
  import { download } from 'matterviz/io'
  import { format_num } from 'matterviz/labels'
  import { ActionMenu, type CmdAction, Icon } from 'svelte-widgets'
  import {
    Code,
    Download,
    Graph,
    Paper,
    PullRequest,
    RSS,
    Unavailable,
  } from 'svelte-widgets/icons'
  import { click_outside, tooltip } from 'svelte-widgets/attachments'
  import type { HTMLAttributes } from 'svelte/elements'
  import {
    ALL_METRICS,
    DISCOVERY_SET_LABELS,
    HYPERPARAMS,
    METADATA_COLS,
  } from '../labels'
  import { assemble_row_data, metric_value } from '../metrics'

  type MetricsRow = ReturnType<typeof assemble_row_data>[number]
  type LinkData = MetricsRow[`Links`]
  type PredFilesDropdown = LinkData[`pred_files`] & { x: number; y: number }
  const export_id = $props.id()
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
    default_sort = DEFAULT_TABLE_SORT,
    sort = $bindable({ ...default_sort }),
    column_preset = `default`,
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    discovery_set?: DiscoverySet
    model_filter?: (model: ModelData) => boolean
    col_filter?: (col: Label) => boolean
    column_labels?: Label[]
    show_row_numbers?: boolean
    filters?: UrlTableFilters
    column_order?: string[]
    default_sort?: SortState
    sort?: SortState
    column_preset?: string
  } = $props()
  let pred_files_dropdown = $state<PredFilesDropdown | null>(null)
  const cohort_models = $derived(ACTIVE_MODELS.filter(model_filter))

  let metrics_data = $derived(
    mark_compared_rows(
      assemble_row_data(discovery_set, model_filter, filters.matches),
      filters.show_selected_only,
    ),
  )
  const mobile = new MediaQuery(`(max-width: 600px)`)
  let mounted = $state(false)
  onMount(() => {
    mounted = true
  })
  const initial_column_order = untrack(() => [...column_order])
  let multi_sort = $derived.by((): { column: string; ascending: boolean }[] => {
    void column_preset
    return []
  })
  const column_defaults = $derived(
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
  const default_column_order = $derived([
    ...new Set([...initial_column_order, ...column_defaults.map((col) => col.id)]),
  ])
  const default_visible_ids = $derived.by(() => {
    const visible = column_defaults.filter((col) => col.visible)
    // Match the server's full table for hydration before applying the mobile view.
    if (!mounted || !mobile.current) return visible.map((col) => col.id)
    const primary = multi_sort[0]?.column ?? sort.column
    const metrics = new Set(
      visible.filter((col) => col.better && col.id !== primary).slice(0, 2),
    )
    return visible
      .filter(
        (col) =>
          col.key === model_name.key ||
          col.id === primary ||
          col.key === model_params.key ||
          metrics.has(col),
      )
      .map((col) => col.id)
  })
  // A new task preset resets custom columns; resizing preserves explicit choices.
  let custom_visible = $derived.by((): string[] | null => {
    void column_preset
    return null
  })
  const columns = $derived(
    column_defaults.map((col) => ({
      ...col,
      default_visible: default_visible_ids.includes(col.id),
      visible: (custom_visible ?? default_visible_ids).includes(col.id),
    })),
  )
  const read_ids = (params: URLSearchParams, key: string): string[] | null => {
    const value = params.get(key)
    if (value === null) return null
    if (value === `none`) return []
    const ids = [...new Set(value.split(`,`))]
    return ids.every((id) => column_defaults.some((col) => col.id === id)) ? ids : null
  }
  bind_url_params(
    (params) => {
      sort = sort_from_query(params, default_sort)
      custom_visible = read_ids(params, `columns`)
      column_order = [
        ...new Set([
          ...(read_ids(params, `column_order`) ?? []),
          ...default_column_order,
        ]),
      ]
      const criteria =
        params
          .get(`multi_sort`)
          ?.split(`,`)
          .filter(Boolean)
          .map((key) => ({
            column: key.startsWith(`-`) ? key.slice(1) : key,
            ascending: !key.startsWith(`-`),
          })) ?? []
      multi_sort = criteria.every(({ column }) =>
        column_defaults.some((col) => col.id === column && col.sortable !== false),
      )
        ? criteria
        : []
    },
    () => [
      ...sort_url_entries(sort, default_sort),
      // Persist mobile defaults as well, so a shared view has the same columns on desktop.
      [
        `columns`,
        custom_visible
          ? custom_visible.join(`,`) || `none`
          : mobile.current
            ? default_visible_ids.join(`,`)
            : ``,
      ],
      [`column_order`, column_order.join(`,`), default_column_order.join(`,`)],
      [
        `multi_sort`,
        multi_sort
          .map(({ column, ascending }) => `${ascending ? `` : `-`}${column}`)
          .join(`,`),
      ],
    ],
  )
  const set_columns = (updated: Column[]): void => {
    const visible = updated.filter((col) => col.visible !== false).map((col) => col.id)
    custom_visible = visible.join(`,`) === default_visible_ids.join(`,`) ? null : visible
  }
  const active_sort = $derived(
    multi_sort.length
      ? multi_sort
      : sort.column
        ? [{ column: sort.column, ascending: sort.dir === `asc` }]
        : [],
  )
  const show_discovery_context = $derived(
    column_labels.some(
      (label) =>
        label.path?.startsWith(`metrics.discovery`) &&
        columns.some((column) => column.key === label.key && column.visible),
    ),
  )
  const show_cps_context = $derived(
    columns.some((column) => column.visible && column.key === `CPS`),
  )
  const cps_total_weight = $derived(
    Object.values(CPS_CONFIG).reduce((total, { weight }) => total + weight, 0),
  )
  async function export_table(export_format: `csv` | `json` | `copy`): Promise<void> {
    const ordered_columns = columns
      .toSorted(
        (left, right) => column_order.indexOf(left.id) - column_order.indexOf(right.id),
      )
      .filter((col) => col.visible)
    const rows = sort_table_rows(
      metrics_data as RowData[],
      active_sort.map(({ column, ascending }) => ({
        key: columns.find((col) => col.id === column)?.key ?? column,
        ascending,
      })),
    ) as MetricsRow[]
    if (export_format !== `json`) {
      const text = table_to_delimited(
        {
          headers: ordered_columns.map(({ label }) => strip_html(label)),
          rows: rows.map((row) =>
            ordered_columns.map(({ id, key = id }) => cell_text(row[key])),
          ),
          numeric: [],
        },
        export_format === `csv` ? `,` : `\t`,
      )
      if (export_format === `copy`) await navigator.clipboard.writeText(text)
      else download(text, `matbench-discovery-${discovery_set}.csv`, `text/csv`)
      return
    }
    const exported_rows = rows.map((row) => {
      const values: Record<string, unknown> = {
        ...row,
        ...Object.fromEntries(
          column_labels.flatMap((label) => {
            const value = metric_value(row.model, label, discovery_set)
            return value === undefined ? [] : [[label.key, value]]
          }),
        ),
        Model: row.model.model_name,
        'Training Set': {
          datasets: row.model.training_sets,
          materials: row.model.n_training_materials,
          structures: row.model.n_training_structures,
        },
        Targets: row.model.targets,
        Org: { logos: row.org_logos, authors: row.authors },
      }
      return Object.fromEntries(
        ordered_columns.map(({ id, key = id }) => {
          const value = values[key]
          return [id, typeof value === `string` ? strip_html(value) : value]
        }),
      )
    })
    const references = Object.fromEntries(
      [
        `wbm_summary`,
        `wbm_initial_atoms`,
        `wbm_relaxed_atoms`,
        `wbm_dft_geo_opt_symprec_1e_2`,
        `wbm_dft_geo_opt_symprec_1e_5`,
        `phonondb_pbe_103_structures`,
        `phonondb_pbe_103_kappa_no_nac`,
        `dynamat_v1_0_md_trajectories`,
        `diatomics_dft_reference`,
      ].map((key) => {
        const entry = data_files[key]
        if (typeof entry === `string`) throw new Error(`Invalid reference file: ${key}`)
        const { url, path, md5 } = entry
        return [key, { url, path, md5 }]
      }),
    )
    download(
      JSON.stringify(
        {
          benchmark_revision: {
            ...BENCHMARK_REVISION,
            development: import.meta.env.DEV,
          },
          url: location.href,
          exported_at: new Date().toISOString(),
          discovery_set,
          cps_discovery_set: `unique_prototypes`,
          filters: { ...filters.as_preset, selected_only: filters.show_selected_only },
          sort: active_sort,
          columns: ordered_columns.map(({ id, key, label, format }) => ({
            id,
            key,
            label,
            format,
          })),
          weights: score_weight_records(),
          references,
          models: rows.map((row) => ({
            model_key: row.model_key,
            model_version: row.model.model_version,
            prediction_files: row.Links.pred_files.files,
          })),
          rows: exported_rows,
        },
        null,
        2,
      ),
      `matbench-discovery-view.json`,
      `application/json`,
    )
  }

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

<div class="ranking-context" aria-label="Ranking context">
  <p>
    <strong>{metrics_data.length} of {cohort_models.length} eligible models</strong>
    · {#if active_sort.length}Sorted by {#each active_sort as criterion, idx}{@const label =
          columns.find(
            (col) => col.id === criterion.column,
          )?.label}{#if idx}{`, then `}{/if}{#if label}{@html label}{:else}{criterion.column}{/if}
        ({criterion.ascending ? `ascending` : `descending`}){/each}{:else}Sort order: see
      column headers{/if}
    {#if show_discovery_context}· Discovery: {DISCOVERY_SET_LABELS[discovery_set]
        .label}{/if}
    {#if show_row_numbers}· Row numbers follow this filtered order.{/if}
  </p>
  {#if show_cps_context}
    <p class="cps-context">
      CPS combines discovery F1 ({format_num(
        CPS_CONFIG.F1.weight / cps_total_weight,
        `.0%`,
      )}), geometry RMSD ({format_num(CPS_CONFIG.RMSD.weight / cps_total_weight, `.0%`)}),
      and phonons κ<sub>SRME</sub> ({format_num(
        CPS_CONFIG.κ_SRME.weight / cps_total_weight,
        `.0%`,
      )}). MD and diatomics are excluded. CPS uses unique-prototype discovery scores.
    </p>
  {/if}
</div>

<HeatmapTable
  data={metrics_data as RowData[]}
  row_key="model_key"
  {columns}
  bind:sort
  bind:multi_sort
  {show_row_numbers}
  default_num_format=".3f"
  bind:show_heatmap={filters.show_heatmap}
  bind:column_order
  on_row_double_click={toggle_row_model}
  {...rest}
  oncontextmenucapture={open_menu}
  class={[`leaderboard`, rest.class]}
  root_style={METRICS_TABLE_ROOT_STYLE}
>
  {#snippet controls()}
    <TableControls
      bind:columns={() => columns, set_columns}
      {filters}
      models={cohort_models}
    >
      {#snippet leading()}
        <a href="/contribute">Submit a model</a>
        <a
          href="/rss.xml"
          title="Follow new model submissions in your RSS reader"
          {@attach tooltip()}
        >
          <Icon icon={RSS} /> RSS
        </a>
      {/snippet}
      {#snippet trailing()}
        <ActionMenu
          aria-labelledby={export_id}
          style="font-size: 0.9em"
          actions={[
            { id: `csv`, label: `CSV`, action: () => export_table(`csv`) },
            {
              id: `json`,
              label: `JSON with provenance`,
              action: () => export_table(`json`),
            },
            { id: `copy`, label: `Copy table`, action: () => export_table(`copy`) },
          ]}
        >
          {#snippet trigger(props)}
            <button
              {...props}
              id={export_id}
              style="margin-left: auto; white-space: nowrap; flex-shrink: 0"
            >
              <Icon icon={Download} /> Export
            </button>
          {/snippet}
        </ActionMenu>
      {/snippet}
    </TableControls>
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
  .ranking-context {
    margin-block: 0.6rem;
    font-size: 0.85em;
    text-align: center;
    color: var(--text-secondary);
    p {
      margin: 0.25em 0;
    }
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
  @media (pointer: coarse), (max-width: 600px) {
    :global(.leaderboard td[data-col='Links']) :is(a, button) {
      display: inline-grid;
      place-items: center;
      min-width: 32px;
      min-height: 36px;
    }
  }
</style>
