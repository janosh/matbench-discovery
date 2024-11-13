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
  export let columns: { label: string; tooltip?: string; style?: string }[] = []
  export let higher_is_better: string[] = []
  export let lower_is_better: string[] = []
  export let sticky_cols: number[] = [0] // default to sticky first column
  export let hide_cols: string[] = [] // just the column labels
  export let format: Record<string, string> = {}
  // set to empty string to hide hint
  export let sort_hint: string = `Click on numerical column headers to sort the table rows by their values`

  const sort_state = writable({ column: ``, ascending: true })

  $: clean_data =
    data?.filter?.((row) => Object.values(row).some((val) => val !== undefined)) ?? []

  function sort_rows(column: string) {
    if ($sort_state.column !== column) {
      $sort_state = {
        column,
        ascending: lower_is_better.includes(column),
      }
    } else {
      $sort_state.ascending = !$sort_state.ascending
    }

    clean_data = clean_data.sort((row1, row2) => {
      const val1 = row1[column]
      const val2 = row2[column]

      if (val1 === val2) return 0
      if (val1 === null || val1 === undefined) return 1
      if (val2 === null || val2 === undefined) return -1

      const modifier = $sort_state.ascending ? 1 : -1
      return val1 < val2 ? -1 * modifier : 1 * modifier
    })
  }

  function calc_color(value: number | string | undefined, col: string) {
    const values = clean_data.map((row) => row[col])
    const range = [min(values) ?? 0, max(values) ?? 1]
    if (lower_is_better.includes(col)) {
      range.reverse()
    }
    const colorScale = scaleSequential()
      .domain(range)
      .interpolator(d3sc.interpolateViridis)

    const bg = colorScale(value)
    const text = choose_bw_for_contrast(null, bg)

    return { bg, text }
  }

  $: visible_columns = columns.filter((col) => !hide_cols.includes(col.label))
</script>

<div class="table-container">
  <table use:titles_as_tooltips>
    <thead>
      <tr>
        {#each visible_columns as { label, tooltip = null, style = null }, col_idx}
          <th on:click={() => sort_rows(label)} title={tooltip} {style}>
            {@html label}
            {#if col_idx == 0 && sort_hint}
              <span title={sort_hint}>
                <Icon icon="octicon:info-16" inline />
              </span>
            {/if}
            {#if $sort_state.column === label}
              <span style="font-size: 0.8em;">
                {$sort_state.ascending ? `↑` : `↓`}
              </span>
            {:else if higher_is_better.includes(label) || lower_is_better.includes(label)}
              <span style="font-size: 0.8em;">
                {higher_is_better.includes(label) ? `↓` : `↑`}
              </span>
            {/if}
          </th>
        {/each}
      </tr>
    </thead>
    <tbody>
      {#each clean_data as row (JSON.stringify(row))}
        <tr animate:flip={{ duration: 500 }}>
          {#each visible_columns as { label, style = null }, col_idx}
            {@const val = row[label]}
            {@const color = calc_color(val, label)}
            <td
              data-col={label}
              data-sort-value={val}
              class:sticky-col={sticky_cols.includes(col_idx)}
              style:background-color={color.bg}
              style:color={color.text}
              {style}
            >
              {#if typeof val === `number` && format[label]}
                {@html pretty_num(val, format[label])}
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
</style>
