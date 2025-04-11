<script lang="ts">
  import { DATASETS, arr_to_str, format_date } from '$lib'
  import HeatmapTable from '$lib/HeatmapTable.svelte'
  import type { Dataset, HeatmapColumn, RowData } from '$lib/types'
  import * as pkg from '$site/package.json'
  import { pretty_num } from 'elementari'

  const columns: HeatmapColumn[] = [
    { label: `Title`, sticky: true },
    {
      label: `Structures`,
      better: `higher`,
      color_scale: `interpolateViridis`,
      scale_type: `log`,
    },
    {
      label: `Materials`,
      better: `higher`,
      color_scale: `interpolateViridis`,
      scale_type: `log`,
    },
    { label: `Created`, tooltip: `Date the dataset was created/started` },
    { label: `Open` },
    { label: `License` },
    { label: `Method` },
    { label: `Links`, sortable: false },
  ]

  // License abbreviations mapped to full names
  const license_map: Record<string, string> = {
    'CC-BY-4.0': `Creative Commons Attribution 4.0 International`,
  }

  const to_spaces = (str: string) => str.replaceAll(`_`, ` `).replaceAll(`-`, ` `)
  const to_title = (str: string) => str.charAt(0).toUpperCase() + str.slice(1)
  const title_case = (str: string) => to_spaces(str).split(` `).map(to_title).join(` `)

  // Process data for table
  const table_data: RowData[] = Object.entries(DATASETS).map(([key, set]) => {
    let { n_structures, n_materials, date_created, license, method, slug } = set
    const license_full = license_map[license]

    const params_tooltip = Object.entries(set.params ?? {})
      .map(([key, value]) => `${title_case(key)}: ${arr_to_str(value)}`)
      .join(`&#013;`)

    const method_str = arr_to_str(method)
    const created_timestamp = date_created ? new Date(date_created).getTime() : null

    return {
      key,
      slug,
      Title: `<a href="/data/${slug}" title="${set.title}">${key}</a>`,
      // Store numeric values directly for coloring, and use cell rendering for display
      Structures: n_structures || null,
      Materials: n_materials || null,
      Open: set.open ? `✅` : `❌`,
      Created: date_created
        ? `<span data-sort-value="${created_timestamp}" title="${format_date(date_created, { weekday: `long` })}">${format_date(date_created)}</span>`
        : `n/a`,
      License: license_full ? `<span title="${license_full}">${license}</span>` : license,
      Method: `<span title="${params_tooltip}">${method_str}</span>`,
      Links: generate_links(set),
    }
  })

  function generate_links(set: Dataset): string {
    const ext_link = (url: string, title: string, text: string) =>
      `<a href="${url}" target="_blank" rel="noopener noreferrer" title="${title}">${text}</a>`
    const links = []

    if (set.url) links.push(ext_link(set.url, `Website`, `🌐`))

    if (set.download_url) links.push(ext_link(set.download_url, `Download`, `💾`))

    if (set.doi) links.push(ext_link(set.doi, `DOI`, `📄`))

    return links.join(` `)
  }

  const yaml_url = `${pkg.repository}/blob/-/data/datasets.yml`
</script>

<svelte:head>
  <title>Training Sets | MatBench Discovery</title>
</svelte:head>

<h1>Datasets</h1>

<p>
  A collection of datasets used for training machine learning models for materials
  discovery.
</p>

<HeatmapTable
  data={table_data}
  {columns}
  initial_sort_column="Created"
  initial_sort_direction="desc"
  fixed_header={true}
  sort_hint=""
>
  {#snippet cell({ col, val })}
    {#if [`Structures`, `Materials`].includes(col.label)}
      {#if val === null}
        n/a
      {:else if typeof val === `number`}
        <span data-sort-value={val} title={val.toLocaleString()}>{pretty_num(val)}</span>
      {/if}
    {:else}
      {@html val}
    {/if}
  {/snippet}
</HeatmapTable>

<p>
  See incorrect data or a dataset that's missing from this list? Suggest an edit to
  <a href={yaml_url} target="_blank" rel="noopener noreferrer">
    {yaml_url.split(`/`).pop()}
  </a>
</p>
