<script lang="ts">
  import { arr_to_str, DATASETS, format_date } from '$lib'
  import { DATASET_METADATA_COLS, title_case } from '$lib/labels'
  import pkg from '$site/package.json'
  import type { RowData } from 'matterviz'
  import { Icon, type IconData } from 'svelte-widgets'
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
  import { HeatmapTable } from 'matterviz'

  const license_map: Record<string, string> = {
    'CC-BY-4.0': `Creative Commons Attribution 4.0 International`,
    'CC-BY-NC-4.0': `Creative Commons Attribution-NonCommercial 4.0 International`,
    MIT: `MIT License`,
  }

  const icon_html = (data: IconData, color?: string): string => {
    const fill = data.fill ?? (data.stroke ? `none` : `currentColor`)
    const stroke = data.stroke ? ` stroke="${data.stroke}"` : ``
    const style = color ? ` style="color:${color}"` : ``
    const inner = `markup` in data ? data.markup : `<path d="${data.d}" />`
    return `<svg fill="${fill}"${stroke}${style} width="1em" height="1em" viewBox="${data.viewBox}">${inner}</svg>`
  }

  const icon_link = (
    href: string | null | undefined,
    title: string,
    glyph: IconData,
    require_http = false,
  ): string | false => {
    if (!href || (require_http && !href.startsWith(`http`))) return false
    return `<a href="${href}" target="_blank" rel="noopener noreferrer" title="${title}" aria-label="${title}">${icon_html(glyph)}</a>`
  }

  const join_links = (links: (string | false)[]): string =>
    links.filter(Boolean).join(` `)

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
    const created_timestamp = date_created ? new Date(date_created).getTime() : null

    const api_links = join_links([
      icon_link(set.native_api, `Native API`, API, true),
      icon_link(set.optimade_api, `OPTIMADE API`, Optimade, true),
    ])
    const resource_links = join_links([
      icon_link(set.url, `Website`, Globe),
      icon_link(set.download_url, `Download`, Download),
      icon_link(set.doi, `DOI`, DOI),
    ])

    return {
      key,
      [DATASET_METADATA_COLS.name.label]:
        `<a href="/data/${slug}" title="${set.name}">${key}</a>`,
      Structures: set.n_structures || null,
      Materials: set.n_materials || null,
      Open: icon_html(
        set.open ? CheckCircle : XCircle,
        set.open ? `lightgreen` : `lightcoral`,
      ),
      Static: icon_html(
        set.static ? CheckCircle : XCircle,
        set.static ? `lightgreen` : `lightcoral`,
      ),
      Created: date_created
        ? `<span data-sort-value="${created_timestamp}" title="${format_date(
            date_created,
            { weekday: `long` },
          )}">${format_date(date_created)}</span>`
        : `n/a`,
      License: license_full ? `<span title="${license_full}">${license}</span>` : license,
      Method: `<span title="${method_str}&#013;${params_tooltip}">${method_str}</span>`,
      API: api_links,
      Links: resource_links,
    }
  })

  const yaml_url = `${pkg.repository}/blob/-/data/datasets.yml`
</script>

<svelte:head>
  <title>Training Sets | MatBench Discovery</title>
</svelte:head>

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
    fixed_header={true}
    sort_hint=""
  />
</section>

<p>
  <Icon icon={Edit} />
  See incorrect data or a dataset that's missing from this list? Suggest an edit to
  <a href={yaml_url} target="_blank" rel="noopener noreferrer">
    {yaml_url.split(`/`).pop()}
  </a>
</p>
