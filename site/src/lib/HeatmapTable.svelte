<script lang="ts">
  import { calc_cell_color } from '$lib/metrics'
  import type { CellSnippetArgs, CellVal, Label, RowData } from '$lib/types'
  import { format_num } from 'matterviz'
  import type { Snippet } from 'svelte'
  import { tooltip } from 'svelte-multiselect/attachments'
  import { flip } from 'svelte/animate'
  import type { HTMLAttributes } from 'svelte/elements'
  import { SvelteMap } from 'svelte/reactivity'

  let {
    data,
    columns = [],
    sort_hint = ``,
    cell,
    special_cells,
    controls,
    initial_sort_column,
    initial_sort_direction,
    fixed_header = false,
    default_num_format = `.3`,
    show_heatmap = $bindable(true),
    heatmap_class = `heatmap`,
    onrowdblclick,
    column_order = $bindable([]),
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    data: RowData[]
    columns?: Label[]
    sort_hint?: string
    cell?: Snippet<[CellSnippetArgs]>
    special_cells?: Record<string, Snippet<[CellSnippetArgs]>>
    controls?: Snippet
    initial_sort_column?: string
    initial_sort_direction?: `asc` | `desc`
    fixed_header?: boolean
    default_num_format?: string
    show_heatmap?: boolean
    heatmap_class?: string
    onrowdblclick?: (event: MouseEvent, row: RowData) => void
    column_order?: string[]
  } = $props()

  // Hacky helper function to detect if a string contains HTML, TODO revisit in future
  function is_html_str(val: unknown): boolean {
    if (typeof val !== `string`) return false
    // Check for common HTML patterns
    return (
      (val.includes(`<`) && val.includes(`>`)) || // Has angle brackets
      val.startsWith(`&lt;`) || // Has HTML entity for <
      val.includes(`<a `) || // Has anchor tag
      val.includes(`<span `) || // Has span tag
      val.includes(`<div `) || // Has div tag
      val.includes(`href=`) || // Has href attribute
      val.includes(`class=`) // Has class attribute
    )
  }

  // Add container reference for binding
  type SortState = { column: string; ascending: boolean }
  let sort_state = $state<SortState>({
    column: initial_sort_column || ``,
    ascending: initial_sort_direction !== `desc`,
  })

  // Helper to make column IDs (needed since column labels in different groups can be repeated)
  const get_col_id = (col: Label) =>
    col.group ? `${col.short ?? col.label} (${col.group})` : (col.short ?? col.label)

  // Initialize and sanitize column_order
  $effect(() => {
    if (columns.length === 0) return
    const ids = columns.map(get_col_id)
    const id_set = new Set(ids)
    if (column_order.length === 0) {
      column_order = ids
    } else {
      const filtered = column_order.filter((id) => id_set.has(id))
      if (filtered.length !== column_order.length || filtered.length !== ids.length) {
        column_order = [...filtered, ...ids.filter((id) => !filtered.includes(id))]
      }
    }
  })

  // Reorder columns based on column_order
  let ordered_columns = $derived.by(() => {
    if (column_order.length === 0) return columns

    const col_map = new SvelteMap(columns.map((col) => [get_col_id(col), col]))
    const ordered: Label[] = []

    // First add columns in the specified order
    for (const col_id of column_order) {
      const col = col_map.get(col_id)
      if (col) {
        ordered.push(col)
        col_map.delete(col_id)
      }
    }

    // Then add any remaining columns that weren't in the order list
    for (const col of col_map.values()) ordered.push(col)

    return ordered
  })

  // Drag state
  let drag_col_id = $state<string | null>(null)
  let drag_over_col_id = $state<string | null>(null)

  function handle_drag_start(event: DragEvent, col: Label) {
    if (!event.dataTransfer) return
    drag_col_id = get_col_id(col)
    event.dataTransfer.effectAllowed = `move`
    event.dataTransfer.setData(`text/html`, ``)
  }

  function handle_drag_over(event: DragEvent, col: Label) {
    event.preventDefault()
    if (!event.dataTransfer) return
    event.dataTransfer.dropEffect = `move`

    // Prevent cross-group drag-over to keep group headers contiguous
    const drag_group = ordered_columns.find((c) => get_col_id(c) === drag_col_id)
      ?.group
    const over_group = col.group
    if (drag_group !== over_group) {
      event.dataTransfer.dropEffect = `none`
      drag_over_col_id = null
      return
    }

    drag_over_col_id = get_col_id(col)
  }

  function handle_drag_leave() {
    drag_over_col_id = null
  }

  function handle_drop(event: DragEvent, target_col: Label) {
    event.preventDefault()

    if (!drag_col_id) return

    // Block cross-group (or group→ungroup) reorders to preserve group contiguity
    const drag_group = ordered_columns.find((c) => get_col_id(c) === drag_col_id)
      ?.group
    if (drag_group !== target_col.group) {
      drag_col_id = null
      drag_over_col_id = null
      return
    }

    const target_col_id = get_col_id(target_col)
    if (drag_col_id === target_col_id) {
      drag_col_id = null
      drag_over_col_id = null
      return
    }

    // Create new order array
    const new_order = [...column_order]
    const drag_idx = new_order.indexOf(drag_col_id)
    const target_idx = new_order.indexOf(target_col_id)

    if (drag_idx === -1 || target_idx === -1) {
      drag_col_id = null
      drag_over_col_id = null
      return
    }

    // Remove dragged column and insert at new position
    new_order.splice(drag_idx, 1)
    new_order.splice(target_idx, 0, drag_col_id)

    column_order = new_order
    drag_col_id = null
    drag_over_col_id = null
  }

  let sorted_data = $derived.by(() => {
    const filtered_data = data?.filter?.((row) =>
      Object.values(row).some((val) => val !== undefined)
    ) ?? []

    if (!sort_state.column) return filtered_data

    const col = ordered_columns.find((c) => get_col_id(c) === sort_state.column)
    if (!col) return filtered_data

    const col_id = get_col_id(col)

    return [...filtered_data].sort((row1, row2) => {
      const val1 = row1[col_id]
      const val2 = row2[col_id]

      if (val1 === val2) return 0

      // Handle null, undefined, and NaN values (always sort to bottom)
      const is_invalid = (val: unknown) =>
        val == null || (typeof val === `number` && Number.isNaN(val))
      if (is_invalid(val1) || is_invalid(val2)) {
        return +is_invalid(val1) - +is_invalid(val2)
      }

      const modifier = sort_state.ascending ? 1 : -1

      // Check if values are HTML strings with data-sort-value attributes
      if (typeof val1 === `string` && typeof val2 === `string`) {
        const sort_val1_match = val1.match(/data-sort-value="([^"]*)"/)
        const sort_val2_match = val2.match(/data-sort-value="([^"]*)"/)

        if (sort_val1_match && sort_val2_match) {
          const sort_val1 = sort_val1_match[1]
          const sort_val2 = sort_val2_match[1]

          // Try to convert to numbers if possible
          const num_val1 = Number(sort_val1)
          const num_val2 = Number(sort_val2)

          if (!isNaN(num_val1) && !isNaN(num_val2)) {
            return num_val1 < num_val2 ? -1 * modifier : 1 * modifier
          }

          // sort strings case-insensitively
          const [lower1, lower2] = [sort_val1.toLowerCase(), sort_val2.toLowerCase()]
          return lower1 > lower2 ? modifier : -1 * modifier
        }
      }

      return (val1 ?? 0) < (val2 ?? 0) ? -1 * modifier : 1 * modifier
    })
  })

  function sort_rows(column: string, group?: string) {
    // Find the column using both label and group if provided
    const col = ordered_columns.find(
      (c) => c.label === column && (c.group === group || c.group === undefined),
    )

    if (!col) return // Skip if column not found
    if (col.sortable === false) return // Skip sorting if column marked as unsortable

    const col_id = get_col_id(col)
    if (sort_state.column !== col_id) {
      sort_state.column = col_id
      sort_state.ascending = col.better === `lower`
    } else sort_state.ascending = !sort_state.ascending
  }

  function calc_color(val: CellVal, col: Label) {
    // Skip color calculation for null values, NaN, or if color_scale is null
    if (
      val === null ||
      val === undefined ||
      col.color_scale === null ||
      typeof val !== `number` ||
      Number.isNaN(val) ||
      !show_heatmap // Disable heatmap colors if show_heatmap is false
    ) return { bg: null, text: null }

    const col_id = get_col_id(col)
    const numeric_vals = sorted_data
      .map((row) => row[col_id])
      .filter((val) => typeof val === `number` && !Number.isNaN(val)) // Type guard to ensure we only get valid numbers

    // Using the shared helper function for color calculation
    return calc_cell_color(
      val,
      numeric_vals,
      col.better === `higher` || col.better === `lower` ? col.better : undefined,
      col.color_scale || `interpolateViridis`,
      col.scale_type || `linear`,
    )
  }

  let visible_columns = $derived(
    ordered_columns.filter((col) => col.visible !== false),
  )

  const sort_indicator = (col: Label, sort_state: SortState) => {
    const col_id = get_col_id(col)
    if (sort_state.column === col_id) {
      // When column is sorted, show ↓ for ascending (smaller values at top)
      // and ↑ for descending (larger values at top)
      return `<span style="font-size: 0.8em;">${
        sort_state.ascending ? `↓` : `↑`
      }</span>`
    } else if (col.better) {
      // When column is not sorted, show arrow indicating which values are better:
      // ↑ for higher-is-better metrics
      // ↓ for lower-is-better metrics
      const sort_dir = col.better === `higher` ? `↑` : `↓`
      return `<span style="font-size: 0.8em;">${sort_dir}</span>`
    }
    return ``
  }
</script>

<div {@attach tooltip()} {...rest} class="table-container {rest.class ?? ``}">
  {#if (sort_state && sort_hint) || controls}
    <div class="table-header">
      {#if sort_state && sort_hint}
        <span class="sort-hint">{sort_hint}</span>
      {/if}
      {#if controls}
        {@render controls()}
      {/if}
    </div>
  {/if}
  <table class:fixed-header={fixed_header} class={heatmap_class} style="grid-column: 2">
    <thead>
      <!-- Don't add a table row for group headers if there are none -->
      {#if visible_columns.some((col) => col.group)}
        <!-- First level headers -->
        <tr class="group-header">
          {#each visible_columns as { label, group, description } (label + group)}
            {#if !group}
              <th></th>
            {:else}
              {@const group_cols = visible_columns.filter((c) => c.group === group)}
              <!-- Only render the group header once for each group by checking if this is the first column of this group -->
              {#if visible_columns.findIndex((c) => c.group === group) ===
            visible_columns.findIndex((c) =>
              c.group === group && c.label === label
            )}
                <th title={description} colspan={group_cols.length}>{@html group}</th>
              {/if}
            {/if}
          {/each}
        </tr>
      {/if}
      <!-- Second level headers -->
      <tr>
        {#each visible_columns as col (col.label + col.group)}
          <th
            title={col.description}
            onclick={() => sort_rows(col.label, col.group)}
            style={col.style}
            class:sticky-col={col.sticky}
            class:not-sortable={col.sortable === false}
            class:dragging={drag_col_id === get_col_id(col)}
            class:drag-over={drag_over_col_id === get_col_id(col)}
            draggable="true"
            aria-dropeffect="move"
            ondragstart={(event: DragEvent & { currentTarget: HTMLElement }) => {
              handle_drag_start(event, col)
              event.currentTarget.setAttribute(`aria-grabbed`, `true`)
            }}
            ondragover={(event) => handle_drag_over(event, col)}
            ondragleave={handle_drag_leave}
            ondrop={(event) => handle_drop(event, col)}
            ondragend={(event: DragEvent & { currentTarget: HTMLElement }) => {
              ;[drag_col_id, drag_over_col_id] = [null, null]
              event.currentTarget.removeAttribute(`aria-grabbed`)
            }}
          >
            {@html col.short ?? col.label}
            {@html sort_indicator(col, sort_state)}
          </th>
        {/each}
      </tr>
    </thead>
    <tbody>
      {#each sorted_data as row (JSON.stringify(row))}
        <tr
          animate:flip={{ duration: 500 }}
          style={row.style}
          class={String(row.class) ?? null}
          ondblclick={onrowdblclick ? (event) => onrowdblclick(event, row) : undefined}
        >
          {#each visible_columns as col (col.label + col.group)}
            {@const val = row[get_col_id(col)]}
            {@const color = calc_color(val, col)}
            <td
              data-col={col.label}
              data-sort-value={is_html_str(val) ? null : val}
              class:sticky-col={col.sticky}
              style:background-color={color.bg}
              style:color={color.text}
              style={col.cell_style ?? col.style}
            >
              {#if special_cells?.[col.label]}
                {@render special_cells[col.label]({ row, col, val })}
              {:else if cell}
                {@render cell({ row, col, val })}
              {:else if typeof val === `number` && !Number.isNaN(val)}
                {format_num(val, col.format ?? default_num_format)}
              {:else if val === undefined || val === null || Number.isNaN(val)}
                <span {@attach tooltip({ content: `Not available` })}>
                  n/a
                </span>
              {:else}
                {@html val}
              {/if}
            </td>
          {/each}
        </tr>
      {/each}
    </tbody>
  </table>
</div>

<style>
  .table-container {
    display: grid;
    grid-template-columns: 1fr min-content 1fr;
    font-size: var(--heatmap-font-size, 0.9em);
  }
  .table-container::-webkit-scrollbar {
    display: none; /* Safari and Chrome */
  }
  table {
    overflow-x: auto;
    overflow-y: hidden;
    /* https://stackoverflow.com/a/38994837 */
    scrollbar-width: none; /* Firefox */
    max-width: 90vw;
  }
  th, td {
    padding: var(--heatmap-cell-padding, 1pt 5pt);
    text-align: var(--heatmap-text-align, left);
    border: var(--heatmap-cell-border, none);
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
  }
  th {
    background: var(--heatmap-header-bg, var(--page-bg));
    position: sticky;
    cursor: pointer;
  }
  th:hover {
    background: var(--heatmap-header-hover-bg, var(--nav-bg));
  }
  th.dragging {
    opacity: 0.4;
    cursor: grabbing;
  }
  th.drag-over {
    outline: 3px solid var(--highlight, #4a9eff);
    outline-offset: -3px;
  }
  th[draggable='true'] {
    cursor: grab;
  }
  .sticky-col {
    position: sticky;
    left: 0;
    background: var(--heatmap-header-bg, var(--page-bg));
    z-index: 1;
  }
  tbody tr:hover {
    filter: var(--heatmap-row-hover-filter, brightness(1.1));
  }
  td[data-sort-value] {
    cursor: default;
  }
  .group-header th {
    border-bottom: 1px solid var(--border);
    text-align: center;
  }
  /* Styles for the table header with sort hint and controls */
  .table-header {
    grid-column: 2;
    width: 100%;
    display: flex;
    place-items: center;
    margin: 10pt auto;
    gap: 2em;
    border-bottom: 1px solid var(--border);
    justify-content: space-between;
  }
  span.sort-hint {
    color: var(--text-muted);
    margin: 0;
  }
  .not-sortable {
    cursor: default;
  }
  tr.highlight {
    background-color: var(--nav-bg) !important;
  }
  tr.highlight, tr.highlight :global(a) {
    color: var(--highlight) !important;
  }
</style>
