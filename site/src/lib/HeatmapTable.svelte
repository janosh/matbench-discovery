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
  import type { HeatmapColumn } from './types'

  type CellVal = string | number | undefined | null
  type RowData = Record<string, CellVal>
  type TableData = RowData[]

  interface Props {
    data: TableData
    columns?: HeatmapColumn[]
    sort_hint?: string
    style?: string | null
    cell?: Snippet<[{ row: RowData; col: HeatmapColumn; val: CellVal }]>
  }

  let {
    data,
    columns = [],
    sort_hint = `Click on column headers to sort table rows`,
    style = null,
    cell,
  }: Props = $props()

  const sort_state = writable({ column: ``, ascending: true })

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
      return val1 < val2 ? -1 * modifier : 1 * modifier
    })
  }

  function calc_color(value: number | string | undefined, col: HeatmapColumn) {
    if (col.color_scale === null || typeof value !== `number`)
      return { bg: null, text: null }

    const col_id = get_col_id(col)
    const values = clean_data.map((row) => row[col_id])
    const range = [min(values) ?? 0, max(values) ?? 1]
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

<div class="table-container" {style}>
  <table use:titles_as_tooltips>
    <thead>
      <!-- Don't add a table row for group headers if there are none -->
      {#if visible_columns.some((col) => col.group)}
        <!-- First level headers -->
        <tr class="group-header">
          {#each visible_columns as col (col.label)}
            {#if !col.group}
              <th></th>
            {:else}
              {@const group_cols = visible_columns.filter((c) => c.group === col.group)}
              {#if columns.indexOf(col) === columns.findIndex((c) => c.group === col.group)}
                <th title={col.tooltip} colspan={group_cols.length}>{@html col.group}</th>
              {/if}
            {/if}
          {/each}
        </tr>
      {/if}
      <!-- Second level headers -->
      <tr>
        {#each visible_columns as col, col_idx (col.label)}
          <th
            title={col.tooltip}
            onclick={() => sort_rows(col.label)}
            style={col.style}
            class:sticky-col={col.sticky}
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
          {#each visible_columns as col (col.label)}
            {@const val = row[get_col_id(col)]}
            {@const color = calc_color(val, col)}
            <td
              data-col={col.label}
              data-sort-value={val}
              class:sticky-col={col.sticky}
              style:background-color={color.bg}
              style:color={color.text}
              style={col.style}
              title={[undefined, null].includes(val) ? `not available` : null}
            >
              {#if cell}
                {@render cell({ row, col, val })}
              {:else if typeof val === `number` && col.format}
                {pretty_num(val, col.format)}
              {:else if [undefined, null].includes(val)}
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
    max-width: 100%;
    scrollbar-width: none;
    margin: auto;
    font-size: var(--heatmap-font-size, 0.9em);
  }

  /* https://stackoverflow.com/a/38994837 */
  .table-container {
    scrollbar-width: none; /* Firefox */
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
</style>
