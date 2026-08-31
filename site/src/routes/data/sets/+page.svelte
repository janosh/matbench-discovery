<script lang="ts">
  import { arr_to_str, DATASETS } from '$lib'
  import { DATASET_METADATA_COLS, title_case } from '$lib/labels'
  import pkg from '$site/package.json'
  import type { CellSnippetArgs, RowData } from 'matterviz'
  import { HeatmapTable } from 'matterviz'
  import { Icon } from 'svelte-widgets'
  import {
    API,
    CheckCircle,
    Databases,
    DOI,
    Download,
    Edit,
    Globe,
    Optimade,
    XCircle,
  } from 'svelte-widgets/icons'

  const license_map: Record<string, string> = {
    'CC-BY-4.0': `Creative Commons Attribution 4.0 International`,
    'CC-BY-NC-4.0': `Creative Commons Attribution-NonCommercial 4.0 International`,
    MIT: `MIT License`,
  }

  // icon links per cell: [key in the cell's link record, tooltip, icon]. Rendered via
  // snippets since HeatmapTable sanitizes HTML strings (dropping inline <svg>)
  const api_links = [
    [`native_api`, `Native API`, API],
    [`optimade_api`, `OPTIMADE API`, Optimade],
  ] as const
  const resource_links = [
    [`url`, `Website`, Globe],
    [`download_url`, `Download`, Download],
    [`doi`, `DOI`, DOI],
  ] as const
  type LinkRecord = Partial<Record<string, string | null>>
  // some API fields hold notes rather than URLs; only link actual http(s) URLs
  const http_url = (url: string | null | undefined): string | null =>
    url?.startsWith(`http`) ? url : null

  const table_data: RowData[] = Object.entries(DATASETS).map(([key, set]) => {
    const { date_created, license, method, slug } = set
    const license_full = license_map[license] ?? license
    const params_tooltip = Object.entries(set.params ?? {})
      .map(
        ([param_key, param_value]) =>
          `${title_case(param_key)}: ${arr_to_str(param_value)}`,
      )
      .join(`&#013;`)
    const method_str = arr_to_str(method)

    return {
      key,
      [DATASET_METADATA_COLS.name.label]:
        `<a href="/data/${slug}" title="${set.name}">${key}</a>`,
      Structures: set.n_structures || null,
      Materials: set.n_materials || null,
      Created: date_created ?? null, // ISO date: HeatmapTable formats and sorts it
      Open: set.open,
      Static: set.static ?? false,
      License: license_full ? `<span title="${license_full}">${license}</span>` : license,
      Method: `<span title="${method_str}&#013;${params_tooltip}">${method_str}</span>`,
      API: {
        native_api: http_url(set.native_api),
        optimade_api: http_url(set.optimade_api),
      },
      Links: { url: set.url, download_url: set.download_url, doi: set.doi },
    }
  })

  const yaml_url = `${pkg.repository}/blob/-/data/datasets.yml`
</script>

<svelte:head>
  <title>Training Sets | MatBench Discovery</title>
</svelte:head>

{#snippet bool_cell({ val }: CellSnippetArgs)}
  <Icon
    icon={val ? CheckCircle : XCircle}
    style="color: {val ? `lightgreen` : `lightcoral`}"
  />
{/snippet}

{#snippet icon_links(
  val: CellSnippetArgs[`val`],
  links: typeof api_links | typeof resource_links,
)}
  {@const record = val as LinkRecord}
  {#each links as [link_key, title, icon] (link_key)}
    {@const href = record[link_key]}
    {#if href}
      <a {href} target="_blank" rel="noopener noreferrer" {title} aria-label={title}>
        <Icon {icon} />
      </a>
    {/if}
  {/each}
{/snippet}

{#snippet api_cell({ val }: CellSnippetArgs)}
  {@render icon_links(val, api_links)}
{/snippet}

{#snippet links_cell({ val }: CellSnippetArgs)}
  {@render icon_links(val, resource_links)}
{/snippet}

<h1>
  <Icon icon={Databases} style="vertical-align: -3pt" /> Datasets
</h1>

<p>
  A collection of datasets used for training machine learning models for materials
  discovery.
</p>
<section class="full-bleed">
  <HeatmapTable
    data={table_data}
    columns={Object.values(DATASET_METADATA_COLS)}
    initial_sort={{ column: `Created`, direction: `desc` }}
    sort_hint=""
    special_cells={{
      [DATASET_METADATA_COLS.open.label]: bool_cell,
      [DATASET_METADATA_COLS.static.label]: bool_cell,
      [DATASET_METADATA_COLS.api.label]: api_cell,
      [DATASET_METADATA_COLS.links.label]: links_cell,
    }}
  />
</section>

<p>
  <Icon icon={Edit} />
  See incorrect data or a dataset that's missing from this list? Suggest an edit to
  <a href={yaml_url} target="_blank" rel="noopener noreferrer">
    {yaml_url.split(`/`).pop()}
  </a>
</p>
