<script lang="ts">
  import { max, min } from 'd3-array'
  import { scaleSequential } from 'd3-scale'
  import * as d3sc from 'd3-scale-chromatic'
  import { choose_bw_for_contrast } from 'elementari/labels'
  import { titles_as_tooltips } from 'svelte-zoo/actions'

  export let data: Record<string, Record<string, unknown>>[] = []
  export let columns: string[] = []
  export let higherIsBetter: Set<string> = new Set()
  export let lowerIsBetter: Set<string> = new Set()

  function cell_bg_color(value: number, col: string): string {
    const values = data.map((d) => d[col])
    const range = [min(values) ?? 0, max(values) ?? 1]
    if (lowerIsBetter.has(col)) {
      range.reverse()
    }
    const colorScale = scaleSequential()
      .domain(range)
      .interpolator(d3sc.interpolateViridis)

    return colorScale(value)
  }

  // Add this type near the top of the file with other types/interfaces
  type ChooseBwActionParams = {
    bgColor?: string | null
    threshold?: number
  }

  // Add this new action before or after the original choose_bw_for_contrast function
  export function choose_bw_for_contrast_action(
    node: HTMLElement,
    params?: ChooseBwActionParams,
  ) {
    const update = (params?: ChooseBwActionParams) => {
      const color = choose_bw_for_contrast(
        node,
        params?.bgColor ?? null,
        params?.threshold,
      )
      node.style.color = color
    }

    update(params) // Initial call

    return { update } // Called when params change
  }
</script>

<div class="table-container">
  <table use:titles_as_tooltips>
    <thead>
      <tr>
        {#each columns as col}
          <th>
            {@html col}
            {#if higherIsBetter.has(col)}
              <span class="arrow">↑</span>
            {:else if lowerIsBetter.has(col)}
              <span class="arrow">↓</span>
            {/if}
          </th>
        {/each}
      </tr>
    </thead>
    <tbody>
      {#each data as row}
        <tr>
          {#each columns as col, idx}
            {@const val = row[col]}
            <td
              data-col={col}
              data-sort-value={val}
              class:sticky-col={idx == 0}
              style:background-color={cell_bg_color(val, col)}
              use:choose_bw_for_contrast_action
              title={val?.toString() ?? ``}
            >
              {@html val?.toString() ?? ``}
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

  .arrow {
    font-size: 0.8em;
    opacity: 0.7;
  }

  tbody tr:hover {
    filter: brightness(1.1);
  }

  td[data-sort-value] {
    cursor: default;
  }
</style>
