<script lang="ts">
  import { max, min } from 'd3-array'
  import { scaleSequential } from 'd3-scale'
  import * as d3sc from 'd3-scale-chromatic'
  import { choose_bw_for_contrast, pretty_num } from 'elementari/labels'
  import { titles_as_tooltips } from 'svelte-zoo/actions'
  import { flip } from 'svelte/animate'
  import { writable } from 'svelte/store'

  export let data: Record<string, Record<string, unknown>>[] = []
  export let columns: string[] = []
  export let higher_is_better: string[] = []
  export let lower_is_better: string[] = []
  export let sep_lines: number[] = []
  export let sticky_cols: number[] = [0] // default to sticky first column
  export let format: Record<string, string> = {}

  const sort_state = writable({ column: ``, ascending: true })
  function sort_rows(column: string) {
    if ($sort_state.column !== column) {
      $sort_state = {
        column,
        ascending: lower_is_better.includes(column),
      }
    } else {
      $sort_state.ascending = !$sort_state.ascending
    }

    data = [...data].sort((row1, row2) => {
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
    const values = data.map((row) => row[col])
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
</script>

<div class="table-container">
  <table use:titles_as_tooltips>
    <thead>
      <tr>
        {#each columns as col}
          <th on:click={() => sort_rows(col)}>
            {@html col}
            {#if higher_is_better.includes(col) || lower_is_better.includes(col)}
              {#if $sort_state.column === col}
                <span style="font-size: 0.8em;">
                  {$sort_state.ascending ? `↑` : `↓`}
                </span>
              {:else}
                <span style="font-size: 0.8em; opacity: 0.5;">
                  {higher_is_better.includes(col) ? `↓` : `↑`}
                </span>
              {/if}
            {:else if $sort_state.column === col}
              <span style="font-size: 0.8em;">
                {$sort_state.ascending ? `↑` : `↓`}
              </span>
            {/if}
          </th>
        {/each}
      </tr>
    </thead>
    <tbody>
      {#each data as row (JSON.stringify(row))}
        <tr animate:flip={{ duration: 500 }}>
          {#each columns as column, idx}
            {@const val = row[column]}
            {@const color = calc_color(val, column)}
            <td
              data-col={column}
              data-sort-value={val}
              class:sticky-col={sticky_cols.includes(idx)}
              class:sep-line={sep_lines.includes(idx)}
              style:background-color={color.bg}
              style:color={color.text}
            >
              {#if typeof val === `number` && format[column]}
                {@html pretty_num(val, format[column])}
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
    width: max-content;
    margin: auto;
  }

  .table-container::-webkit-scrollbar {
    display: none;
  }

  table {
    border-collapse: collapse;
    width: 100%;
    font-size: 0.9rem;
  }

  th,
  td {
    padding: 1pt 3pt;
    text-align: left;
    border: none;
    white-space: nowrap;
    max-width: 300px;
    overflow: hidden;
    text-overflow: ellipsis;
  }

  th {
    background: var(--night);
    position: sticky;
    top: 0;
    z-index: 1;
    cursor: pointer;
    user-select: none;
  }

  th:hover {
    background: var(--night-lighter, #2a2a2a);
  }

  .sep-line {
    border-right: 1px solid black;
  }

  .sticky-col {
    position: sticky;
    left: 0;
    background: var(--night);
    z-index: 2;
  }
  tr:nth-child(odd) td.sticky-col {
    background: rgb(15, 14, 14);
  }

  tbody tr:hover {
    filter: brightness(1.1);
  }

  td[data-sort-value] {
    cursor: default;
  }
</style>
