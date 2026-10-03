<script lang="ts">
  import { goto } from '$app/navigation'
  import data_files from '$pkg/data-files.yml'
  import { bind_url_params, type UrlTableFilters } from '#lib/url-state.svelte.js'
  import { sort_from_query, sort_url_entries } from 'matterviz/url-params'
  import { MediaQuery } from 'svelte/reactivity'
  import { onMount, untrack } from 'svelte'
  import { CPS_CONFIG } from '#lib/combined-scores.svelte.js'
  import OrgLogos from '#lib/model/OrgLogos.svelte'
  import TableControls from '#lib/table/TableControls.svelte'
  import { append_better_hint, missing_metric_reason } from '#lib/metrics.js'
  import {
    comparison,
    mark_compared_rows,
    row_model_key,
    toggle_row_model,
  } from '#lib/model-comparison.svelte.js'
  import {
    ACTIVE_MODELS,
    make_table_filters,
    MODELS,
    score_weight_records,
  } from '#lib/models.svelte.js'
  import type { DiscoverySet, Label, ModelData } from '#lib/types.js'
  import type {
    CellSnippet,
    CellSnippetArgs,
    Column,
    RowData,
    TableSort,
  } from 'matterviz/table'
  import {
    cell_text,
    HeatmapTable,
    is_invalid,
    table_to_delimited,
  } from 'matterviz/table'
  import { escape_html, strip_html } from 'matterviz/utils'
  import { download } from 'matterviz/io'
  import { format_num } from 'matterviz/labels'
  import { ActionMenu, type CmdAction, Icon, Popover } from 'svelte-widgets'
  import {
    Code,
    Download,
    Graph,
    Paper,
    PullRequest,
    RSS,
    Unavailable,
  } from 'svelte-widgets/icons'
  import { tooltip } from 'svelte-widgets/attachments'
  import type { HTMLAttributes } from 'svelte/elements'
  import {
    ALL_METRICS,
    DISCOVERY_SET_LABELS,
    HYPERPARAMS,
    METADATA_COLS,
    plain_label,
  } from '../labels.js'
  import { assemble_row_data, metric_value } from '../metrics.js'

  type MetricsRow = ReturnType<typeof assemble_row_data>[number]
  type LinkData = MetricsRow[`Links`]
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
    default_sort = { column: `CPS`, dir: `desc` }, // omitted from shared URLs
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
    default_sort?: TableSort
    sort?: TableSort
    column_preset?: string
  } = $props()
  // rows left after search/filters in the table's displayed order (all pages), bound
  // from HeatmapTable so export and comparison seeding follow what the table shows
  let visible_rows = $state<RowData[]>([])
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
    column_labels.map((col): Column => ({
      ...col,
      id: col.group ? `${col.key} (${col.group})` : col.key,
      cell: cells[col.key],
      ...(column_labels === default_columns && {
        color_scale: heatmap_disabled_cols.has(col.key) ? null : col.color_scale,
        ...(col === model_name && { style: `padding-left: 0;${col.style ?? ``}` }),
      }),
      better: col.better ?? undefined, // null (no direction) isn't a Column value
      description: append_better_hint(col),
      visible: col.visible !== false && col_filter(col),
    })),
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
  async function export_table(
    export_format: `csv` | `json` | `copy` | `copy_html`,
  ): Promise<void> {
    const ordered_columns = columns
      .toSorted(
        (left, right) => column_order.indexOf(left.id) - column_order.indexOf(right.id),
      )
      .filter((col) => col.visible)
    const rows = visible_rows as MetricsRow[]
    if (export_format !== `json`) {
      const matrix = {
        headers: ordered_columns.map(({ label }) => plain_label(label)),
        rows: rows.map((row) =>
          ordered_columns.map(({ id, key = id }) => cell_text(row[key])),
        ),
        numeric: [],
      }
      const text = table_to_delimited(matrix, export_format === `csv` ? `,` : `\t`)
      if (export_format === `copy_html`) {
        const html_row = (values: string[], tag: `th` | `td`) =>
          `<tr>${values.map((value) => `<${tag} style="border: 1px solid #ccc; padding: 4px 8px; text-align: left">${escape_html(value)}</${tag}>`).join(``)}</tr>`
        const html = `<table style="border-collapse: collapse"><thead>${html_row(matrix.headers, `th`)}</thead><tbody>${matrix.rows.map((row) => html_row(row, `td`)).join(``)}</tbody></table>`
        await navigator.clipboard.write([
          new ClipboardItem({
            'text/html': new Blob([html], { type: `text/html` }),
            'text/plain': new Blob([text], { type: `text/plain` }),
          }),
        ])
      } else if (export_format === `copy`) await navigator.clipboard.writeText(text)
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
          filters: { ...filters.config, selected_only: filters.show_selected_only },
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

{#snippet affiliation_cell({ row }: CellSnippetArgs)}
  {@const metrics_row = row as MetricsRow}
  <OrgLogos org_logos={metrics_row.org_logos} authors={metrics_row.authors} />
{/snippet}

{#snippet metric_cell({ row, col, val }: CellSnippetArgs)}
  {@const coverage =
    col.key === ALL_METRICS.pbe_vib_freq_error.key
      ? (row as MetricsRow).model.metrics?.diatomics?.pbe_vib_freq_coverage
      : undefined}
  <!-- the fit count only earns its width when some elements lack a valid fit -->
  {#if coverage && coverage.n_valid < coverage.n_eligible}
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
  <Popover class="pred-files-dropdown" aria-label="Files for {links.pred_files.name}">
    {#snippet trigger(trigger_props)}
      <button
        style="background: none; padding: 0"
        aria-label="Download model prediction files"
        {...trigger_props}
      >
        <Icon icon={Graph} />
      </button>
    {/snippet}
    <h4 style="margin: 0">Files for {links.pred_files.name}</h4>
    <ol style="margin: 0; padding-left: 1em">
      {#each links.pred_files.files as { name: file_name, url } (url)}
        <li>
          <a href={url} target="_blank" rel="noopener noreferrer">{@html file_name}</a>
        </li>
      {/each}
    </ol>
  </Popover>
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
  bind:visible_rows
  {show_row_numbers}
  default_num_format=".3f"
  bind:show_heatmap={filters.show_heatmap}
  bind:column_order
  on_row_double_click={toggle_row_model}
  {...rest}
  oncontextmenucapture={open_menu}
  class={[`leaderboard`, rest.class]}
  root_style="--heatmap-sticky-cell-odd-bg: linear-gradient(var(--table-odd), var(--table-odd)), var(--page-bg); --heatmap-row-num-padding-left: 0; --heatmap-column-max-width: 14.4em;"
>
  {#snippet controls()}
    <TableControls
      bind:columns={() => columns, set_columns}
      {filters}
      models={cohort_models}
      {visible_rows}
    >
      {#snippet leading()}
        <a href="/contribute">Submit a model</a>
        <a
          href="/rss.xml"
          data-toolbar-optional
          style="margin-inline: 0.25em"
          title="Follow new model submissions in your RSS reader"
          {@attach tooltip()}
        >
          <Icon icon={RSS} /> RSS
        </a>
        <button
          class="table-guide"
          data-toolbar-optional
          type="button"
          title={`Select a column heading to sort. Hover labels for definitions and n/a cells for missing-result explanations. Click plotted models or double-click table rows to highlight them across tables and plots. Use Compare to view them side by side.\n\nTraining Set counts distinct materials, with relaxation frames in parentheses. When only frame counts are available, those are shown instead.`}
          {@attach tooltip({
            touch_focus: true,
            placement: `bottom`,
            wrap: `normal`,
            style: `white-space: pre-line; --tooltip-max-width: 32rem; --tooltip-padding: 0.75em 1em`,
          })}>How to read the table</button
        >
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
            {
              id: `copy_html`,
              label: `Copy table as HTML`,
              action: () => export_table(`copy_html`),
            },
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

<style>
  .table-guide {
    font: inherit;
    background: none;
    color: var(--link-color);
    padding: 0;
    cursor: help;
  }
  .ranking-context {
    margin-block: 0.6rem;
    font-size: 0.85em;
    text-align: center;
    color: var(--text-secondary);
    p {
      margin: 0.25em 0;
    }
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
