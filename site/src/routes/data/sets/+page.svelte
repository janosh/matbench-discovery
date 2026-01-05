<script lang="ts">
  import { arr_to_str, DATASETS, format_date, Icon } from '$lib'
  import { icon_data } from '$lib/icons'
  import { DATASET_METADATA_COLS, title_case } from '$lib/labels'
  import type { RowData } from '$lib/types'
  import pkg from '$site/package.json'
  import { HeatmapTable } from 'matterviz'

  const license_map: Record<string, string> = {
    'CC-BY-4.0': `Creative Commons Attribution 4.0 International`,
    'CC-BY-NC-4.0': `Creative Commons Attribution-NonCommercial 4.0 International`,
    MIT: `MIT License`,
  }

  const icon = (name: keyof typeof icon_data, color?: string) => {
    const data = icon_data[name]
    const fill = `stroke` in data ? `none` : `currentColor`
    const stroke = `stroke` in data ? `stroke="currentColor"` : ``
    return `<svg fill="${fill}" ${stroke} ${
      color ? `style="color:${color}"` : ``
    } width="1em" height="1em" viewBox="${data.viewBox}">
    ${data.path.trim().startsWith(`<`) ? data.path : `<path d="${data.path}" />`}
    </svg>`
  }

  const table_data: RowData[] = Object.entries(DATASETS).map(([key, set]) => {
    const { date_created, license, method, slug } = set
    const license_full = license_map[license] ?? license
    const params_tooltip = Object.entries(set.params ?? {}).map(([k, v]) =>
      `${title_case(k)}: ${arr_to_str(v)}`
    ).join(`&#013;`)
    const method_str = arr_to_str(method)
    const created_timestamp = date_created ? new Date(date_created).getTime() : null

    const api_links = [
      set.native_api?.startsWith(`http`) &&
      `<a href="${set.native_api}" target="_blank" rel="noopener noreferrer" title="Native API">${
        icon(`API`)
      }</a>`,
      set.optimade_api?.startsWith(`http`) &&
      `<a href="${set.optimade_api}" target="_blank" rel="noopener noreferrer" title="OPTIMADE API">${
        icon(`Optimade`)
      }</a>`,
    ].filter(Boolean).join(` `)

    const resource_links = [
      set.url &&
      `<a href="${set.url}" target="_blank" title="Website">${icon(`Globe`)}</a>`,
      set.download_url &&
      `<a href="${set.download_url}" target="_blank" title="Download">${
        icon(`Download`)
      }</a>`,
      set.doi &&
      `<a href="${set.doi}" target="_blank" title="DOI">${icon(`DOI`)}</a>`,
    ].filter(Boolean).join(` `)

    return {
      key,
      [DATASET_METADATA_COLS.name.label]:
        `<a href="/data/${slug}" title="${set.name}">${key}</a>`,
      Structures: set.n_structures || null,
      Materials: set.n_materials || null,
      Open: icon(
        set.open ? `CheckCircle` : `XCircle`,
        set.open ? `lightgreen` : `lightcoral`,
      ),
      Static: icon(
        set.static ? `CheckCircle` : `XCircle`,
        set.static ? `lightgreen` : `lightcoral`,
      ),
      Created: date_created
        ? `<span data-sort-value="${created_timestamp}" title="${
          format_date(date_created, { weekday: `long` })
        }">${format_date(date_created)}</span>`
        : `n/a`,
      License: license_full
        ? `<span title="${license_full}">${license}</span>`
        : license,
      Method:
        `<span title="${method_str}&#013;${params_tooltip}">${method_str}</span>`,
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
  <Icon icon="Databases" style="vertical-align: -3pt" /> Datasets
</h1>

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
  <Icon icon="Edit" />
  See incorrect data or a dataset that's missing from this list? Suggest an edit to
  <a href={yaml_url} target="_blank" rel="noopener noreferrer">
    {yaml_url.split(`/`).pop()}
  </a>
</p>
