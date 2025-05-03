<script lang="ts">
  import { DATASETS, arr_to_str, format_date } from '$lib'
  import HeatmapTable from '$lib/HeatmapTable.svelte'
  import { DATASET_METADATA_COLS, title_case } from '$lib/labels'
  import type { Dataset, RowData } from '$lib/types'
  import pkg from '$site/package.json'

  // License abbreviations mapped to full names
  const license_map: Record<string, string> = {
    'CC-BY-4.0': `Creative Commons Attribution 4.0 International`,
    'CC-BY-NC-4.0': `Creative Commons Attribution-NonCommercial 4.0 International`,
    MIT: `MIT License`,
  }

  const ext_link = (url: string, title: string, text: string) =>
    `<a href="${url}" target="_blank" rel="noopener noreferrer" title="${title}">${text}</a>`

  // Process data for table
  const table_data: RowData[] = Object.entries(DATASETS).map(([key, set]) => {
    let { date_created, license, method, slug } = set
    const license_full = license_map[license] ?? license

    const params_tooltip = Object.entries(set.params ?? {})
      .map(([key, value]) => `${title_case(key)}: ${arr_to_str(value)}`)
      .join(`&#013;`) // Use HTML entity for newline in title attribute

    const method_str = arr_to_str(method)
    const created_timestamp = date_created ? new Date(date_created).getTime() : null

    return {
      key,
      slug,
      [DATASET_METADATA_COLS.name.label]:
        `<a href="/data/${slug}" title="${set.name}">${key}</a>`,
      // Store numeric values directly for coloring, and use cell rendering for display
      Structures: set.n_structures || null,
      Materials: set.n_materials || null,
      Open: `<svg color="${set.open ? `lightgreen` : `lightcoral`}"><use href="#icon-${set.open ? `check` : `x`}" /></svg>`,
      Static: `<svg color="${set.static ? `lightgreen` : `lightcoral`}"><use href="#icon-${set.static ? `check` : `x`}" /></svg>`,
      Created: date_created
        ? `<span data-sort-value="${created_timestamp}" title="${format_date(date_created, { weekday: `long` })}">${format_date(date_created)}</span>`
        : `n/a`,
      License: license_full ? `<span title="${license_full}">${license}</span>` : license,
      Method: `<span title="${method_str}&#013;${params_tooltip}">${method_str}</span>`,
      API: generate_api_links(set),
      Links: generate_resource_links(set),
    }
  })

  function generate_api_links(set: Dataset): string {
    let links = ``
    if (set.native_api?.startsWith(`http`))
      links +=
        `&nbsp;` +
        ext_link(set.native_api, `Native API`, `<svg><use href="#icon-api" /></svg>`)

    if (set.optimade_api?.startsWith(`http`))
      links +=
        `&nbsp;` +
        ext_link(
          set.optimade_api,
          `OPTIMADE API`,
          `<svg><use href="#icon-optimade" /></svg>`,
        )

    return links
  }

  function generate_resource_links(set: Dataset): string {
    let links = ``
    if (set.url)
      links +=
        `&nbsp;` + ext_link(set.url, `Website`, `<svg><use href="#icon-globe" /></svg>`)
    if (set.download_url)
      links +=
        `&nbsp;` +
        ext_link(set.download_url, `Download`, `<svg><use href="#icon-download" /></svg>`)
    if (set.doi)
      links += `&nbsp;` + ext_link(set.doi, `DOI`, `<svg><use href="#icon-doi" /></svg>`)
    return links
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
<section class="full-bleed">
  <HeatmapTable
    data={table_data}
    columns={Object.values(DATASET_METADATA_COLS)}
    initial_sort_column="Created"
    initial_sort_direction="desc"
    fixed_header={true}
    sort_hint=""
  />
</section>

<p>
  See incorrect data or a dataset that's missing from this list? Suggest an edit to
  <a href={yaml_url} target="_blank" rel="noopener noreferrer">
    {yaml_url.split(`/`).pop()}
  </a>
</p>
