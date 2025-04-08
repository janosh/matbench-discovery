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
    { label: `Links`, sortable: false },
  ]

  type TableRow = {
    key: string
    Title: string
    Structures: string
    Materials: string
    Open: string
    Created: string
    License: string
    Method: string
    Links: string
    sort_values: Record<string, string | number | null>
    [key: string]: string | Record<string, string | number | null>
  }

  // License abbreviations mapped to full names
  const license_map: Record<string, string> = {
    'CC-BY-4.0': `Creative Commons Attribution 4.0 International`,
  }

  const span_wrap = (val: string | number) =>
    `<span data-sort-value="${val}" title="${val.toLocaleString()}">${pretty_num(Number(val))}</span>`
  const to_spaces = (str: string) => str.replaceAll(`_`, ` `).replaceAll(`-`, ` `)
  const to_title = (str: string) => str.charAt(0).toUpperCase() + str.slice(1)
  const title_case = (str: string) => to_spaces(str).split(` `).map(to_title).join(` `)

  // convert array types to strings and handle missing values
  function arr_to_str(method: unknown): string {
    if (!method) return `n/a`
    if (Array.isArray(method)) return method.join(`, `)
    return String(method)
  }

  // Process data for table
  const datasets = Object.entries(DATASETS).map(([key, set]) => {
    let { n_structures, n_materials, date_created, license } = set
    const license_full = license_map[license]

    const params_tooltip = Object.entries(set.params ?? {})
      .map(([key, value]) => `${title_case(key)}: ${arr_to_str(value)}`)
      .join(`&#013;`)

    const method_str = arr_to_str(set.params?.method)
    const created_timestamp = date_created ? new Date(date_created).getTime() : null

    return {
      key,
      Title: `<a href="${set.url}" title="${set.title}">${key}</a>`,
      Structures: n_structures ? span_wrap(n_structures) : `n/a`,
      Materials: n_materials ? span_wrap(n_materials) : `n/a`,
      Open: set.open ? `✅` : `❌`,
      Created: date_created
        ? `<span data-sort-value="${created_timestamp}" title="${format_date(date_created, { weekday: `long` })}">${format_date(date_created)}</span>`
        : `n/a`,
      License: license_full ? `<span title="${license_full}">${license}</span>` : license,
      Method: `<span title="${params_tooltip}">${method_str}</span>`,
      Links: generate_links(set),
      sort_values: {
        Title: key,
        Structures: n_structures ?? null,
        Materials: n_materials ?? null,
        Open: set.open ? 1 : 0,
        Created: created_timestamp,
        License: license,
        Method: method_str,
      },
    }
  }) as TableRow[]

  function generate_links(set: Dataset): string {
    const ext_link = (url: string, title: string, text: string) =>
      `<a href="${url}" target="_blank" rel="noopener noreferrer" title="${title}">${text}</a>`
    const links = []

    if (set.url) links.push(ext_link(set.url, `Website`, `🌐`))

    if (set.download_url) links.push(ext_link(set.download_url, `Download`, `💾`))

    if (set.doi) links.push(ext_link(set.doi, `DOI`, `📄`))

    return links.join(` `)
  }

  let sort_state = {
    column: `Created`,
    ascending: true,
  }

  // Initial sort
  let table_data = sort_rows(`Created`)

  function sort_rows(column: string): TableRow[] {
    const col = columns.find((c) => c.label === column)
    if (!col) return datasets

    // Skip sorting if column is explicitly marked as not sortable
    if (col.sortable === false) return datasets

    if (sort_state.column !== column) {
      sort_state = {
        column,
        ascending: col.better === `lower`,
      }
    } else {
      sort_state.ascending = !sort_state.ascending
    }

    return datasets.sort((row1, row2) => {
      const val1 = row1.sort_values[column]
      const val2 = row2.sort_values[column]

      if (val1 === val2) return 0
      if (val1 === null || val1 === undefined) return 1
      if (val2 === null || val2 === undefined) return -1

      const modifier = sort_state.ascending ? 1 : -1

      return val1 < val2 ? -1 * modifier : 1 * modifier
    })
  }

  const yaml_url = `https://github.com/janosh/matbench-discovery/blob/main/data/datasets.yml`
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

See a dataset that's missing from this list or incorrect data? Suggest an edit to
<a href={yaml_url} target="_blank" rel="noopener noreferrer"
  >{yaml_url.split(`/`).pop()}</a
>

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
