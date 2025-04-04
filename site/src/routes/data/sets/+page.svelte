<script lang="ts">
  import { DATASETS } from '$lib'
  import { format_date } from '$lib/metrics-table-helpers'
  import type { Dataset, HeatmapColumn } from '$lib/types'
  import { pretty_num } from 'elementari'
  import 'iconify-icon'
  import { titles_as_tooltips } from 'svelte-zoo/actions'
  import { flip } from 'svelte/animate'

  const columns: HeatmapColumn[] = [
    { label: `Title`, sticky: true },
    { label: `Structures`, better: `higher` },
    { label: `Materials`, better: `higher` },
    { label: `Open` },
    { label: `Created`, better: `higher` },
    { label: `License` },
    { label: `Method` },
    { label: `Code` },
    { label: `Links`, sortable: false },
  ]

  type TableRow = {
    Title: string
    Structures: string
    Materials: string
    Open: string
    Created: string
    License: string
    Method: string
    Code: string
    Links: string
    [key: string]: string // Index signature to allow string indexing
  }

  const span_wrap = (val: string | number) =>
    `<span data-sort-value="${val}" title="${val.toLocaleString()}">${pretty_num(Number(val))}</span>`

  // Process data for table
  const datasets = Object.entries(DATASETS).map(([key, set]) => {
    let { n_structures, n_materials, date_created } = set
    return {
      key,
      Title: `<span title="${set.title}">${key}</span>`,
      Structures: n_structures ? span_wrap(n_structures) : `n/a`,
      Materials: n_materials ? span_wrap(n_materials) : `n/a`,
      Open: set.open ? `✅` : `❌`,
      Created: date_created
        ? `<span data-sort-value="${new Date(date_created).getTime()}">${format_date(date_created)}</span>`
        : `n/a`,
      License: set.license || `n/a`,
      Method: set.params?.method || `n/a`,
      Code: set.params?.code || `n/a`,
      Links: generate_links(set),
    }
  })

  function generate_links(set_data: Dataset): string {
    const ext_link = (url: string, title: string, text: string) =>
      `<a href="${url}" target="_blank" rel="noopener noreferrer" title="${title}">${text}</a>`
    const links = []

    if (set_data.url) links.push(ext_link(set_data.url, `Website`, `🌐`))

    if (set_data.download_url)
      links.push(ext_link(set_data.download_url, `Download`, `💾`))

    if (set_data.doi) links.push(ext_link(set_data.doi, `DOI`, `📄`))

    return links.join(` `)
  }

  let sort_state = {
    column: `Created`,
    ascending: true,
  }

  // Initial sort
  let table_data: TableRow[] = sort_rows(`Created`)

  function sort_rows(column: string): TableRow[] {
    const col = columns.find((c) => c.label === column)
    if (!col) return [...datasets]

    // Skip sorting if column is explicitly marked as not sortable
    if (col.sortable === false) return [...datasets]

    if (sort_state.column !== column) {
      sort_state = {
        column,
        ascending: col.better === `lower`,
      }
    } else {
      sort_state.ascending = !sort_state.ascending
    }

    return [...datasets].sort((row1, row2) => {
      const val1 = row1[column]
      const val2 = row2[column]

      if (val1 === val2) return 0
      if (val1 === null || val1 === undefined || val1 === `n/a`) return 1
      if (val2 === null || val2 === undefined || val2 === `n/a`) return -1

      const modifier = sort_state.ascending ? 1 : -1

      // Check if values are HTML strings with data-sort-value attributes
      if (typeof val1 === `string` && typeof val2 === `string`) {
        const match1 = val1.match(/data-sort-value="([^"]*)"/)
        const match2 = val2.match(/data-sort-value="([^"]*)"/)

        if (match1 && match2) {
          const sort_val1 = match1[1]
          const sort_val2 = match2[1]

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

  // Generate sort indicator
  function sort_indicator(column: string): string {
    const col = columns.find((c) => c.label === column)
    if (!col) return ``

    if (sort_state.column === column) {
      // When column is sorted, show ↓ for ascending (smaller values at top)
      // and ↑ for descending (larger values at top)
      return `<span style="font-size: 0.8em;">${sort_state.ascending ? `↓` : `↑`}</span>`
    } else if (col.better) {
      // When column is not sorted, show arrow indicating which values are better:
      // ↑ for higher-is-better metrics
      // ↓ for lower-is-better metrics
      return `<span style="font-size: 0.8em;">${col.better === `higher` ? `↑` : `↓`}</span>`
    }
    return ``
  }
</script>

<svelte:head>
  <title>Training Sets | MatBench Discovery</title>
</svelte:head>

<h1>Datasets</h1>

A collection of datasets used for training machine learning models for materials
discovery.

<div class="table-wrapper">
  <div class="table-container" use:titles_as_tooltips>
    <table class="training-sets-table">
      <thead>
        <tr>
          {#each columns as col (col.label)}
            <th
              onclick={() => {
                if (col.sortable != false) table_data = sort_rows(col.label)
              }}
              class:sticky-col={col.sticky}
              class:sortable={col.sortable != false}
            >
              {@html col.label}
              {@html sort_indicator(col.label)}
              {#if col.label === `Title`}
                <span title="Click on column headers to sort">
                  <iconify-icon icon="octicon:info-16" inline></iconify-icon>
                </span>
              {/if}
            </th>
          {/each}
        </tr>
      </thead>
      <tbody>
        {#each table_data as row (row.key)}
          <tr animate:flip={{ duration: 500 }}>
            {#each columns as col (col.label)}
              <td class:sticky-col={col.sticky} data-col={col.label}>
                {@html row[col.label]}
              </td>
            {/each}
          </tr>
        {/each}
      </tbody>
    </table>
  </div>
</div>

<style>
  /* Use negative margin technique for full width, matching +page.svelte */
  .table-wrapper {
    width: calc(100vw - 40px);
    margin: 1em calc(-50vw + 50% + 20px);
    display: flex;
    justify-content: center;
    overflow-x: auto;
  }
  table {
    border-collapse: collapse;
  }
  :is(th, td) {
    padding: 3pt 10pt;
    text-align: left;
    border: none;
  }
  tbody tr:nth-child(odd) {
    background: rgba(0, 0, 0, 0.5);
  }
  tbody tr:hover {
    filter: brightness(1.1);
  }
  th.sortable {
    cursor: pointer;
  }
</style>
