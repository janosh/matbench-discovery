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
    // Array of column IDs to control display order. IDs are derived as:
    // - Ungrouped columns: col.short ?? col.label
    // - Grouped columns: `${col.short ?? col.label} (${col.group})`
    // This allows persisting/restoring column order across sessions.
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

    // Initialize if empty, otherwise sync: keep valid IDs in order, append new ones
    const filtered = column_order.filter((id) => id_set.has(id))
    const needs_update = column_order.length === 0 ||
      filtered.length !== column_order.length ||
      filtered.length !== ids.length

    // Guard against unnecessary churn: skip update if order is already identical
    if (
      needs_update &&
      !(column_order.length && filtered.length === ids.length &&
        column_order.every((val, idx) => val === ids[idx]))
    ) {
      column_order = column_order.length === 0
        ? ids
        : [...filtered, ...ids.filter((id) => !filtered.includes(id))]
    }
  })

  // Reorder columns based on column_order
  let ordered_columns = $derived.by(() => {
    if (column_order.length === 0) return columns

    const col_map = new SvelteMap(columns.map((col) => [get_col_id(col), col]))

    // Add columns in specified order, then any remaining columns that weren't in the order list
    const ordered = column_order
      .map((id) => col_map.get(id))
      .filter(Boolean) as Label[]

    const ordered_ids = new Set(ordered.map(get_col_id))
    const remaining = columns.filter((col) => !ordered_ids.has(get_col_id(col)))

    return [...ordered, ...remaining]
  })

  let drag_col_id = $state<string | null>(null)
  let drag_over_col_id = $state<string | null>(null)

  // Returns 'left' or 'right' to indicate which side of target to insert dragged column
  function get_drag_side(target_col_id: string) {
    if (!drag_col_id) return null
    const drag_idx = column_order.indexOf(drag_col_id)
    const target_idx = column_order.indexOf(target_col_id)
    if (drag_idx === -1 || target_idx === -1) return null
    return drag_idx < target_idx ? `right` : `left`
  }

  function reset_drag_state() {
    drag_col_id = null
    drag_over_col_id = null
  }

  const get_drag_col_group = () =>
    ordered_columns.find((col) => get_col_id(col) === drag_col_id)?.group

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
    if (get_drag_col_group() !== col.group) {
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

    // Block cross-group (or group→ungroup) reorders to preserve group contiguity
    if (!drag_col_id || drag_col_id === get_col_id(target_col)) {
      reset_drag_state()
      return
    }

    // Block cross-group reorders to preserve group contiguity
    if (get_drag_col_group() !== target_col.group) {
      reset_drag_state()
      return
    }

    const target_col_id = get_col_id(target_col)
    const drag_idx = column_order.indexOf(drag_col_id)
    const target_idx = column_order.indexOf(target_col_id)

    if (drag_idx === -1 || target_idx === -1) {
      reset_drag_state()
      return
    }

    // Reorder: remove dragged column, then insert at target position
    // After removal, target_idx naturally points to the correct insertion spot
    const new_order = [...column_order]
    new_order.splice(drag_idx, 1)
    new_order.splice(target_idx, 0, drag_col_id)
    column_order = new_order
    reset_drag_state()
  }

  let sorted_data = $derived.by(() => {
    const filtered_data = data?.filter?.((row) =>
      Object.values(row).some((val) => val !== undefined)
    ) ?? []

    if (!sort_state.column) return filtered_data

    const col = ordered_columns.find((col) => get_col_id(col) === sort_state.column)
    if (!col) return filtered_data

    const col_id = get_col_id(col)
    const modifier = sort_state.ascending ? 1 : -1

    // Helper to check if value is invalid (null, undefined, NaN)
    const is_invalid = (val: unknown) =>
      val == null || (typeof val === `number` && Number.isNaN(val))

    return [...filtered_data].sort((row1, row2) => {
      const val1 = row1[col_id]
      const val2 = row2[col_id]

      if (val1 === val2) return 0

      // Push invalid values to bottom
      if (is_invalid(val1) || is_invalid(val2)) {
        return +is_invalid(val1) - +is_invalid(val2)
      }

      // Handle HTML strings with data-sort-value attributes
      if (typeof val1 === `string` && typeof val2 === `string`) {
        const match1 = val1.match(/data-sort-value="([^"]*)"/)
        const match2 = val2.match(/data-sort-value="([^"]*)"/)

        if (match1 && match2) {
          const [sort_val1, sort_val2] = [match1[1], match2[1]]
          const [num1, num2] = [Number(sort_val1), Number(sort_val2)]

          // Use numeric comparison if both parse as numbers
          if (!isNaN(num1) && !isNaN(num2)) {
            if (num1 === num2) return 0
            return num1 < num2 ? -modifier : modifier
          }

          // Otherwise sort strings using localeCompare for natural ordering
          return sort_val1.localeCompare(sort_val2, undefined, {
            numeric: true,
            sensitivity: `base`,
          }) * modifier
        }
      }

      return (val1 ?? 0) < (val2 ?? 0) ? -modifier : modifier
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
    const is_sorted = sort_state.column === col_id

    // Show ↓ for ascending/↑ for descending when sorted
    // Show ↑ for higher-is-better/↓ for lower-is-better when not sorted
    const arrow = is_sorted
      ? (sort_state.ascending ? `↓` : `↑`)
      : (col.better === `higher` ? `↑` : col.better === `lower` ? `↓` : ``)

    return arrow ? `<span style="font-size: 0.8em;">${arrow}</span>` : ``
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
          {@const col_id = get_col_id(col)}
          {@const drag_side = drag_over_col_id === col_id
            ? get_drag_side(col_id)
            : null}
          <th
            title={col.description}
            onclick={() => {
              if (!drag_col_id) sort_rows(col.label, col.group)
            }}
            style={col.style}
            class:sticky-col={col.sticky}
            class:not-sortable={col.sortable === false}
            class:dragging={drag_col_id === col_id}
            data-drag-side={drag_side}
            draggable="true"
            aria-dropeffect="move"
            aria-sort={sort_state.column === col_id
            ? (sort_state.ascending ? `ascending` : `descending`)
            : `none`}
            ondragstart={(event: DragEvent & { currentTarget: HTMLElement }) => {
              handle_drag_start(event, col)
              event.currentTarget.setAttribute(`aria-grabbed`, `true`)
            }}
            ondragover={(event) => handle_drag_over(event, col)}
            ondragleave={handle_drag_leave}
            ondrop={(event) => handle_drop(event, col)}
            ondragend={(event: DragEvent & { currentTarget: HTMLElement }) => {
              reset_drag_state()
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
  th[data-drag-side='left'] {
    border-left: 4px solid var(--highlight, #4a9eff);
  }
  th[data-drag-side='right'] {
    border-right: 4px solid var(--highlight, #4a9eff);
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
