<script lang="ts">
  import DATASETS from '$data/datasets.yml'
  import { arr_to_str } from '$lib'
  import { title_case } from '$lib/labels'
  import { ACTIVE_MODELS } from '$lib/models.svelte'
  import DynamicScatter from '$lib/plot/DynamicScatter.svelte'
  import {
    bind_url_params,
    sort_from_query,
    sort_url_entries,
    type SortState,
  } from '$lib/url-state.svelte'
  import pkg from '$site/package.json'
  import type { CellSnippetArgs, Column } from 'matterviz/table'
  import { HeatmapTable, middle_ellipsis_parts } from 'matterviz/table'
  import { format_num } from 'matterviz/labels'
  import { Icon, Popover } from 'svelte-widgets'
  import {
    API,
    Databases,
    DOI,
    Download,
    Edit,
    Globe,
    Info,
    Optimade,
  } from 'svelte-widgets/icons'
  import { valid_query_param } from 'svelte-widgets/url-params'

  const colors = {
    access: { public: `#25836d`, partial: `#b16c00`, unreleased: `#9771bb` },
    role: {
      training: `#397ec5`,
      validation: `#b16c00`,
      test: `#25836d`,
      repository: `#9771bb`,
    },
  }
  const access_descriptions = {
    public: `The corpus counted here is publicly available; license conditions still apply.`,
    partial: `Only a subset is public. Open the dataset details for release coverage.`,
    unreleased: `The training corpus is not publicly released.`,
  }
  const display_value = (value: unknown): string =>
    typeof value === `string` ? value : arr_to_str(value)
  const api_links = [
    [`native_api`, `Native API`, API],
    [`optimade_api`, `OPTIMADE API`, Optimade],
  ] as const
  const resource_links = [
    [`url`, `Website`, Globe],
    [`download_url`, `Download`, Download],
    [`doi`, `DOI`, DOI],
  ] as const
  const datasets = Object.entries(DATASETS).map(([key, dataset]) => ({
    ...dataset,
    key,
    n_models: ACTIVE_MODELS.filter((model) => model.training_sets.includes(key)).length,
    release:
      dataset.static == null ? `Unknown` : dataset.static ? `Fixed release` : `Updated`,
    method_label:
      [dataset.method, dataset.params?.functional]
        .filter(Boolean)
        .map(display_value)
        .join(` · `) || `Unknown`,
  }))
  type DatasetRow = (typeof datasets)[number]
  const column_defs = {
    Name: {
      key: `key`,
      cell: name_cell,
      description: `Name of the dataset`,
      sticky: true,
    },
    Structures: {
      key: `n_structures`,
      cell: size_cell,
      description: `Number of structures in the dataset. Any system with atomic positions and energy/force/stress labels is counted as a structure incl. successive ionic steps in MD/geometry optimization trajectories.`,
      better: `higher`,
      scale_type: `log`,
      format: `.3s`,
    },
    Materials: {
      key: `n_materials`,
      cell: size_cell,
      description: `Distinct materials or molecules. Dataset notes explain the counting scope.`,
      better: `higher`,
      scale_type: `log`,
      format: `.3s`,
    },
    Created: {
      key: `date_created`,
      description: `Date the dataset was created/started`,
    },
    Access: {
      key: `access`,
      cell: category_cell,
      description: `Availability of the full corpus counted here. Public access does not imply unrestricted licensing.`,
    },
    Role: {
      key: `role`,
      cell: category_cell,
      description: `Primary purpose of this entry; subsets may be used for other purposes.`,
    },
    Models: {
      key: `n_models`,
      cell: models_cell,
      color_scale: null,
      description: `Active models explicitly declaring this training set. Composite datasets count under their own names, without inferring use of every source structure.`,
    },
    Static: {
      key: `release`,
      label: `Release`,
      description: `Fixed release or continuously updated repository. Dates mark release or repository creation, not the last update.`,
      style: `text-align: center;`,
    },
    License: {
      key: `license`,
      description: `License under which the dataset is published`,
    },
    Method: {
      key: `method_label`,
      cell: method_cell,
      description: `Method and exchange-correlation functional. Details list the code and other calculation settings.`,
      style: `max-width: 12em;`,
    },
    API: {
      cell: links_cell,
      description: `API docs (OPTIMADE or native)`,
      sortable: false,
    },
    Links: {
      cell: links_cell,
      description: `Relevant links for the dataset`,
      sortable: false,
    },
  } satisfies Record<string, Partial<Column<DatasetRow>>>
  const columns: Column<DatasetRow>[] = Object.entries(column_defs).map(
    ([id, column]) => ({ id, label: id, ...column }),
  )
  const filter_options = [
    [`access`, new Set(Object.keys(colors.access))],
    [`role`, new Set(Object.keys(colors.role))],
    [`method`, new Set(datasets.flatMap(({ method }) => method ?? []).toSorted())],
  ] as const
  const sortable_columns = new Set(
    columns.filter((col) => col.sortable !== false).map(({ id }) => id),
  )
  const default_sort = { column: `Models`, dir: `desc` as const }
  let sort = $state<SortState>({ ...default_sort })
  const default_filters = { q: ``, access: ``, role: ``, method: `` }
  let filters = $state({ ...default_filters })
  const plot_options = (
    [`Created`, `Structures`, `Materials`, `Models`, `Access`, `Role`] as const
  ).map((label) => ({
    key: column_defs[label].key,
    label,
    description: column_defs[label].description,
    ...(label === `Access` || label === `Role`
      ? { categories: colors[column_defs[label].key] }
      : { format: label === `Created` ? `%Y` : label === `Models` ? `d` : `.3~s` }),
  }))
  const default_plot = {
    x: `date_created`,
    y: `n_structures`,
    color: `access`,
    size: `n_models`,
  }
  let plot = $state({ ...default_plot })
  const plot_dims = [`x`, `y`, `color`, `size`] as const
  const plot_keys = new Set(plot_options.map(({ key }) => key))
  const numeric_keys = new Set(
    plot_options.filter((option) => !(`categories` in option)).map(({ key }) => key),
  )
  const read_url = (params: URLSearchParams) => {
    filters.q = params.get(`q`) ?? ``
    for (const [key, options] of filter_options) {
      filters[key] = valid_query_param(params, key, ``, options)
    }
    for (const dim of plot_dims) {
      plot[dim] = valid_query_param(
        params,
        dim,
        default_plot[dim],
        dim === `color` ? plot_keys : numeric_keys,
      )
    }
    sort = sort_from_query(params, default_sort, sortable_columns)
  }
  bind_url_params(read_url, () => [
    ...Object.entries(filters),
    ...plot_dims.map((dim) => [dim, plot[dim], default_plot[dim]] as const),
    ...sort_url_entries(sort, default_sort),
  ])
  const search_words = $derived(filters.q.toLowerCase().trim().split(/\s+/))
  const filtered = $derived(
    datasets.filter((dataset) => {
      const search_text =
        `${dataset.key} ${dataset.name} ${dataset.description}`.toLowerCase()
      return (
        filter_options.every(
          ([key]) =>
            !filters[key] ||
            [dataset[key]].flat().some((value) => value === filters[key]),
        ) && search_words.every((word) => search_text.includes(word))
      )
    }),
  )
  const has_filters = $derived(Object.values(filters).some(Boolean))
  const reset_filters = () => {
    filters = { ...default_filters }
  }
</script>

<svelte:head>
  <title>Datasets | Matbench Discovery</title>
</svelte:head>

{#snippet name_cell({ row: dataset }: CellSnippetArgs<DatasetRow>)}
  <a href="/data/{dataset.slug}" title={dataset.name}>{dataset.key}</a>
  <Popover
    class="dataset-info"
    trigger_mode="hover"
    trap_focus={false}
    aria-labelledby="about-{dataset.slug}"
    style="max-width: min(34rem, 90vw); white-space: normal; font-weight: normal"
  >
    {#snippet trigger(trigger_props)}
      <button
        class="info"
        type="button"
        aria-labelledby="about-{dataset.slug}"
        {...trigger_props}
        onclick={trigger_props.onmouseenter}
        ><Icon icon={Info} /><span id="about-{dataset.slug}" hidden
          >About {dataset.key}</span
        ></button
      >
    {/snippet}
    <strong>{dataset.name}</strong>
    {@html dataset.description_html ?? ``}
    {#if dataset.contains?.length}
      <p>
        Derived from: {#each dataset.contains as source, idx (source)}{#if idx > 0},
          {/if}<a href="/data/{DATASETS[source].slug}">{source}</a>{/each}
      </p>
    {/if}
    {#each Object.entries(dataset.notes_html ?? {}) as [title, note] (title)}
      <strong>{title}</strong>
      {@html note}
    {/each}
    <a href="/data/{dataset.slug}">Dataset details →</a>
  </Popover>
{/snippet}

{#snippet size_cell({ val, row, col }: CellSnippetArgs<DatasetRow>)}
  <span
    title={val == null
      ? col.id === `Structures`
        ? `Structure count not reported; see ${row.key} details for coverage.`
        : `Unique material count not reported; structures may include multiple frames of the same material.`
      : Number(val).toLocaleString()}
    >{val == null ? `n/a` : format_num(Number(val), `.3~s`)}</span
  >
{/snippet}

{#snippet category_cell({ row, col, val }: CellSnippetArgs<DatasetRow>)}
  <span title={col.id === `Access` ? access_descriptions[row.access] : undefined}>
    {title_case(String(val))}
  </span>
{/snippet}

{#snippet method_cell({ row }: CellSnippetArgs<DatasetRow>)}
  {@const [start, end] = middle_ellipsis_parts(row.method_label)}
  <span
    class="method"
    title={`${row.method_label}\n${
      Object.entries(row.params ?? {})
        .map(([key, value]) => `${title_case(key)}: ${display_value(value)}`)
        .join(`\n`) || `Calculation settings not reported.`
    }`}><span>{start}</span><span>{end}</span></span
  >
{/snippet}

{#snippet models_cell({ row }: CellSnippetArgs<DatasetRow>)}
  {#if row.n_models}
    <a
      href="/?{new URLSearchParams({ train: row.key, targets: `` })}"
      aria-label="View {row.n_models} models trained on {row.key}">{row.n_models}</a
    >
  {:else}
    <span title="No active models explicitly declare this training set.">0</span>
  {/if}
{/snippet}

{#snippet links_cell({ row, col }: CellSnippetArgs<DatasetRow>)}
  <span style="display: inline-flex; gap: 0.25em">
    {#each col.id === `API` ? api_links : resource_links as [link_key, title, icon] (link_key)}
      {@const href = row[link_key]}
      {#if href}
        <a
          {href}
          target="_blank"
          rel="noopener noreferrer"
          {title}
          aria-label="{title} for {row.key}"><Icon {icon} /></a
        >
      {/if}
    {/each}
  </span>
{/snippet}

<h1 id="datasets"><Icon icon={Databases} style="vertical-align: -3pt" /> Datasets</h1>
<p>
  Explore materials datasets for training, validation, and testing, plus source
  repositories. Open a dataset’s info button for coverage and caveats.
</p>

<section class="full-bleed">
  <div class="filters" aria-label="Dataset filters" style="margin-bottom: 1.5rem">
    <input
      type="search"
      aria-label="Search datasets"
      placeholder="Search datasets…"
      bind:value={filters.q}
    />
    {#each filter_options as [key, options] (key)}
      <label>
        {title_case(key)}
        <select aria-label={title_case(key)} bind:value={filters[key]}>
          <option value="">All</option>
          {#each options as value (value)}
            <option {value}>{key === `method` ? value : title_case(value)}</option>
          {/each}
        </select>
      </label>
    {/each}
    <span role="status">{filtered.length} of {datasets.length} datasets</span>
    {#if has_filters}<button type="button" onclick={reset_filters}>Reset filters</button
      >{/if}
  </div>
  <HeatmapTable
    data={filtered}
    {columns}
    bind:sort
    initial_sort={{ column: default_sort.column, direction: default_sort.dir }}
    row_key="key"
    row_animation_ms={300}
    sort_hint=""
  />
  {#if filtered.length === 0}<p class="empty">
      No datasets match these filters. <button type="button" onclick={reset_filters}
        >Show all datasets</button
      >
    </p>{/if}
</section>

<section id="dataset-growth">
  <h2 style="text-align: center">
    {plot.x === default_plot.x && plot.y === default_plot.y
      ? `Dataset Sizes Over Time`
      : `${plot_options.find(({ key }) => key === plot.y)?.label} vs ${plot_options.find(({ key }) => key === plot.x)?.label}`}
  </h2>
  <p>
    Sizes refer to the listed corpus, including unreleased portions. Repository dates mark
    creation; their sizes are later snapshots. Entries missing a selected value are
    omitted. Click a point to open its dataset.
  </p>
  <DynamicScatter
    models={filtered}
    item_name="datasets"
    get_identity={({ key, slug }) => ({ key, name: key, href: `/data/${slug}` })}
    options={plot_options}
    bind:x_key={plot.x}
    bind:y_key={plot.y}
    bind:color_key={plot.color}
    bind:size_key={plot.size}
    legend={null}
    style="height: 460px"
  />
</section>

<p>
  <Icon icon={Edit} /> See incorrect data or a missing dataset? Suggest an edit to
  <a
    href="{pkg.repository}/blob/main/data/datasets.yml"
    target="_blank"
    rel="noopener noreferrer">datasets.yml</a
  >.
</p>

<style>
  input,
  select,
  .filters button {
    font: inherit;
    padding: 0.2em 0.4em;
    border: 1px solid var(--border);
    border-radius: 4px;
    background: var(--page-bg);
    color: inherit;
  }
  input::placeholder {
    color: inherit;
    opacity: 1;
  }
  .filters {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    justify-content: center;
    font-size: 0.85rem;
    color: var(--text-color);
    gap: 0.4em 0.75em;
    input {
      flex: 1;
      min-width: 12em;
      max-width: 22em;
    }
    label {
      display: flex;
      align-items: center;
      gap: 0.4em;
    }
  }
  .method {
    display: flex;
    span:first-child {
      overflow: hidden;
      text-overflow: ellipsis;
    }
    span:last-child {
      flex-shrink: 0;
    }
  }
  button.info {
    background: none;
    border: none;
    padding: 0.2em;
    margin-left: 0.25em;
    color: var(--link-color);
  }
  :global(.dataset-info :is(p, ul, ol)) {
    margin-block: 0.35em;
  }
  .empty {
    text-align: center;
  }
</style>
