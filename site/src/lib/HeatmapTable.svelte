<script lang="ts">
  import Icon from '@iconify/svelte'
  import { max, min } from 'd3-array'
  import { scaleSequential } from 'd3-scale'
  import * as d3sc from 'd3-scale-chromatic'
  import { choose_bw_for_contrast, pretty_num } from 'elementari/labels'
  import { titles_as_tooltips } from 'svelte-zoo/actions'
  import { flip } from 'svelte/animate'
  import { writable } from 'svelte/store'

  type TableData = Record<string, string | number | undefined>[]

  export let data: TableData
  export let columns: {
    group?: string
    label: string
    tooltip?: string
    style?: string
  }[] = []
  export let higher_is_better: string[] = []
  export let lower_is_better: string[] = []
  export let sticky_cols: number[] = [0] // default to sticky first column
  export let hide_cols: string[] = [] // just the column labels
  export let format: Record<string, string> = {}
  // set to empty string to hide hint
  export let sort_hint: string = `Click on column headers to sort table rows`
  export let style: string | null = null

  const sort_state = writable({ column: ``, ascending: true })

  $: clean_data =
    data?.filter?.((row) => Object.values(row).some((val) => val !== undefined)) ?? []

  // Helper to make column IDs (needed since column labels in different groups can be repeated)
  const get_col_id = (col: { group?: string; label: string }) =>
    col.group ? `${col.label} (${col.group})` : col.label

  function sort_rows(column: string) {
    const col = columns.find((c) => c.label === column)
    const col_id = get_col_id(col)

    if ($sort_state.column !== col_id) {
      $sort_state = {
        column: col_id,
        ascending: lower_is_better.includes(col_id),
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

  function calc_color(
    value: number | string | undefined,
    col: { group?: string; label: string },
  ) {
    const col_id = get_col_id(col)
    const values = clean_data.map((row) => row[col_id])
    const range = [min(values) ?? 0, max(values) ?? 1]
    if (lower_is_better.includes(col_id)) {
      range.reverse()
    }

    const color_scale = scaleSequential()
      .domain(range)
      .interpolator(d3sc.interpolateViridis)

    const bg = color_scale(value)
    const text = choose_bw_for_contrast(null, bg)

    return { bg, text }
  }

  $: visible_columns = columns.filter((col) => !hide_cols.includes(col.label))

  const sort_indicator = (col: { group?: string; label: string }) => {
    const col_id = get_col_id(col)
    if ($sort_state.column === col_id) {
      return `<span style="font-size: 0.8em;">${$sort_state.ascending ? `↑` : `↓`}</span>`
    } else if (higher_is_better.includes(col_id) || lower_is_better.includes(col_id)) {
      return `<span style="font-size: 0.8em;">${
        higher_is_better.includes(col_id) ? `↑` : `↓`
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
          {#each visible_columns as col}
            {#if !col.group}
              <th />
            {:else}
              {@const group_cols = visible_columns.filter((c) => c.group === col.group)}
              {#if columns.indexOf(col) === columns.findIndex((c) => c.group === col.group)}
                <th colspan={group_cols.length}>{@html col.group}</th>
              {/if}
            {/if}
          {/each}
        </tr>
      {/if}
      <!-- Second level headers -->
      <tr>
        {#each visible_columns as col, col_idx}
          <th on:click={() => sort_rows(col.label)} style={col.style}>
            {@html col.label}
            {@html sort_indicator(col)}
            {#if col_idx == 0 && sort_hint}
              <span title={sort_hint}>
                <Icon icon="octicon:info-16" inline />
              </span>
            {/if}
          </th>
        {/each}
      </tr>
    </thead>
    <tbody>
      {#each clean_data as row (JSON.stringify(row))}
        <tr animate:flip={{ duration: 500 }}>
          {#each visible_columns as col, col_idx}
            {@const val = row[get_col_id(col)]}
            {@const color = calc_color(val, col)}
            <td
              data-col={col.label}
              data-sort-value={val}
              class:sticky-col={sticky_cols.includes(col_idx)}
              style:background-color={color.bg}
              style:color={color.text}
              style={col.style}
            >
              {#if typeof val === `number` && format[get_col_id(col)]}
                {@html pretty_num(val, format[get_col_id(col)])}
              {:else if [undefined, null].includes(val)}
                <span title="not available">n/a</span>
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
    z-index: 2;
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
