<script lang="ts">
  import { max, min } from 'd3-array'
  import { scaleSequential } from 'd3-scale'
  import * as d3sc from 'd3-scale-chromatic'
  import { choose_bw_for_contrast, pretty_num } from 'elementari/labels'
  import 'iconify-icon'
  import type { Snippet } from 'svelte'
  import { titles_as_tooltips } from 'svelte-zoo/actions'
  import { flip } from 'svelte/animate'
  import { writable } from 'svelte/store'
  import type { CellVal, HeatmapColumn, RowData, TableData } from './types'

  interface Props {
    data: TableData
    columns?: HeatmapColumn[]
    sort_hint?: string
    style?: string | null
    cell?: Snippet<[{ row: RowData; col: HeatmapColumn; val: CellVal }]>
    controls?: Snippet
    initial_sort_column?: string
    initial_sort_direction?: `asc` | `desc`
    fixed_header?: boolean
  }

  let {
    data,
    columns = [],
    sort_hint = `Click on column headers to sort table rows`,
    style = null,
    cell,
    controls,
    initial_sort_column,
    initial_sort_direction,
    fixed_header = false,
  }: Props = $props()

  // Add container reference for binding
  let container: HTMLDivElement

  const sort_state = writable({
    column: initial_sort_column || ``,
    ascending: initial_sort_direction !== `desc`,
  })

  let clean_data = $state(data)
  $effect(() => {
    clean_data =
      data?.filter?.((row) => Object.values(row).some((val) => val !== undefined)) ?? []
  })

  // Helper to make column IDs (needed since column labels in different groups can be repeated)
  const get_col_id = (col: HeatmapColumn) =>
    col.group ? `${col.label} (${col.group})` : col.label

  function sort_rows(column: string) {
    const col = columns.find((c) => c.label === column)
    if (!col) return // Skip if column not found

    // Skip sorting if column is explicitly marked as not sortable
    if (col.sortable === false) return

    const col_id = get_col_id(col)

    if ($sort_state.column !== col_id) {
      $sort_state = {
        column: col_id,
        ascending: col.better === `lower`,
      }
    } else {
      $sort_state.ascending = !$sort_state.ascending
    }

    clean_data = clean_data.sort((row1, row2) => {
      const val1 = row1[col_id]
      const val2 = row2[col_id]

      if (val1 === val2) return 0
      if (val1 === null || val1 === undefined) return 1
      if (val2 === null || val2 === undefined) return -1

      const modifier = $sort_state.ascending ? 1 : -1

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

          return sort_val1 < sort_val2 ? -1 * modifier : 1 * modifier
        }
      }

      return val1 < val2 ? -1 * modifier : 1 * modifier
    })
  }

  function calc_color(value: CellVal, col: HeatmapColumn) {
    // Skip color calculation for null values or if color_scale is null
    if (
      value === null ||
      value === undefined ||
      col.color_scale === null ||
      typeof value !== `number`
    ) {
      return { bg: null, text: null }
    }

    const col_id = get_col_id(col)
    const numericValues = clean_data
      .map((row) => row[col_id])
      .filter((val): val is number => typeof val === `number`) // Type guard to ensure we only get numbers

    if (numericValues.length === 0) {
      return { bg: null, text: null }
    }

    const range = [min(numericValues) ?? 0, max(numericValues) ?? 1]
    if (col.better === `lower`) {
      range.reverse()
    }

    // Use custom color scale if specified, otherwise fall back to viridis
    const scale_name = col.color_scale || `interpolateViridis`
    const interpolator = d3sc[scale_name] || d3sc.interpolateViridis

    const color_scale = scaleSequential().domain(range).interpolator(interpolator)

    const bg = color_scale(value)
    const text = choose_bw_for_contrast(null, bg)

    return { bg, text }
  }

  let visible_columns = $derived(columns.filter((col) => !col.hidden))

  const sort_indicator = (
    col: HeatmapColumn,
    sort_state: { column: string; ascending: boolean },
  ) => {
    const col_id = get_col_id(col)
    if (sort_state.column === col_id) {
      // When column is sorted, show ↓ for ascending (smaller values at top)
      // and ↑ for descending (larger values at top)
      return `<span style="font-size: 0.8em;">${sort_state.ascending ? `↓` : `↑`}</span>`
    } else if (col.better) {
      // When column is not sorted, show arrow indicating which values are better:
      // ↑ for higher-is-better metrics
      // ↓ for lower-is-better metrics
      return `<span style="font-size: 0.8em;">${
        col.better === `higher` ? `↑` : `↓`
      }</span>`
    }
    return ``
  }
</script>

<!-- Table header with sort hint and controls side by side -->
<div class="table-header">
  {#if Object.keys($sort_state).length && sort_hint}
    <div class="sort-hint">{sort_hint}</div>
  {/if}

  <!-- Add controls rendering here -->
  {#if controls}
    <div class="controls-container">
      {@render controls()}
    </div>
  {/if}
</div>

<div bind:this={container} class="table-container" {style}>
  <table use:titles_as_tooltips class:fixed-header={fixed_header} class="heatmap">
    <thead>
      <!-- Don't add a table row for group headers if there are none -->
      {#if visible_columns.some((col) => col.group)}
        <!-- First level headers -->
        <tr class="group-header">
          {#each visible_columns as { label, group, tooltip } (label + group)}
            {#if !group}
              <th></th>
            {:else}
              {@const group_cols = visible_columns.filter((c) => c.group === group)}
              <!-- Only render the group header once for each group by checking if this is the first column of this group -->
              {#if columns.findIndex((c) => c.label === label) === columns.findIndex((c) => c.group === group)}
                <th title={tooltip} colspan={group_cols.length}>{@html group}</th>
              {/if}
            {/if}
          {/each}
        </tr>
      {/if}
      <!-- Second level headers -->
      <tr>
        {#each visible_columns as col, col_idx (col.label + col.group)}
          <th
            title={col.tooltip}
            onclick={() => sort_rows(col.label)}
            style={col.style}
            class:sticky-col={col.sticky}
            class:not-sortable={col.sortable === false}
          >
            {@html col.label}
            {@html sort_indicator(col, $sort_state)}
            {#if col_idx == 0 && sort_hint}
              <span title={sort_hint}>
                <iconify-icon icon="octicon:info-16" inline></iconify-icon>
              </span>
            {/if}
          </th>
        {/each}
      </tr>
    </thead>
    <tbody>
      {#each clean_data as row (JSON.stringify(row))}
        <tr animate:flip={{ duration: 500 }}>
          {#each visible_columns as col (col.label + col.group)}
            {@const val = row[get_col_id(col)]}
            {@const color = calc_color(val, col)}
            <td
              data-col={col.label}
              data-sort-value={val}
              class:sticky-col={col.sticky}
              style:background-color={color.bg}
              style:color={color.text}
              style={col.style}
              title={typeof val === `undefined` || val === null ? `not available` : null}
            >
              {#if cell}
                {@render cell({ row, col, val })}
              {:else if typeof val === `number` && col.format}
                {pretty_num(val, col.format)}
              {:else if val === undefined || val === null}
                n/a
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
    overflow-x: auto;
    margin: auto;
    font-size: var(--heatmap-font-size, 0.9em);
    /* https://stackoverflow.com/a/38994837 */
    scrollbar-width: none; /* Firefox */
    max-width: 90vw;
  }
  .table-container::-webkit-scrollbar {
    display: none; /* Safari and Chrome */
  }

  th,
  td {
    padding: var(--heatmap-cell-padding, 1pt 3pt);
    text-align: var(--heatmap-text-align, left);
    border: var(--heatmap-cell-border, none);
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
  }

  th {
    background: var(--heatmap-header-bg, var(--night));
    position: sticky;
    cursor: pointer;
  }

  th:hover {
    background: var(--heatmap-header-hover-bg, var(--night-lighter, #2a2a2a));
  }

  .sticky-col {
    position: sticky;
    left: 0;
    background: var(--heatmap-header-bg, var(--night));
    z-index: 1;
  }
  tr:nth-child(odd) td.sticky-col {
    background: var(--heatmap-row-odd-bg, rgb(15, 14, 14));
  }

  tbody tr:hover {
    filter: var(--heatmap-row-hover-filter, brightness(1.1));
  }

  td[data-sort-value] {
    cursor: default;
  }

  .group-header th {
    border-bottom: 1px solid black;
    text-align: center;
  }

  /* Styles for the table header with sort hint and controls */
  .table-header {
    display: flex;
    align-items: center;
    justify-content: space-between;
    margin-bottom: 0.5rem;
    flex-wrap: wrap;
    gap: 0.5rem;
    padding: 0.25rem 0;
    border-bottom: 1px solid rgba(255, 255, 255, 0.1);
  }

  .sort-hint {
    font-size: 0.85em;
    color: var(--text-muted, #aaa);
    margin: 0;
  }

  .controls-container {
    display: inline-flex;
    align-items: center;
    margin-left: auto;
  }

  .not-sortable {
    cursor: default;
  }
</style>
