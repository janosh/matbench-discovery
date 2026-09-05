<script lang="ts">
  import { arr_to_str, DATASETS, format_date } from '$lib'
  import { format_relative_time, title_case } from '$lib/labels'
  import { format_num } from 'matterviz'
  import { Icon } from 'svelte-widgets'
  import {
    Calendar,
    CalendarPlus,
    Code,
    DOI,
    Database,
    Download,
    Email,
    GitHub,
    Globe,
    Lattice,
    License,
    Lock,
    ORCID,
    Unlock,
  } from 'svelte-widgets/icons'
  import pkg from '$site/package.json'
  import { tooltip } from 'svelte-widgets/attachments'
  import type { PageData } from './$types'
  import MPtrjTargetDistros from './MPtrjTargetDistros.svelte'

  let { data }: { data: PageData } = $props()
  let dataset = $derived(data.dataset)
  const link_props = { target: `_blank`, rel: `noopener noreferrer` }

  let created_ago = $derived(format_relative_time(dataset.date_created))
  let added_ago = $derived(format_relative_time(dataset.date_added))
  let dataset_links = $derived([
    [dataset.url, `Website`, Globe, `View dataset website`],
    [dataset.download_url, `Download`, Download, `Download dataset`],
    [dataset.doi, `DOI`, DOI, `Digital Object Identifier`],
    [
      `${pkg.repository}/blob/main/data/datasets.yml`,
      `Source`,
      Code,
      `View source YAML file`,
    ],
  ] as const)

  // params object -> [label, value] pairs
  const format_params = (
    params: Record<string, unknown> | undefined,
  ): [string, string][] =>
    Object.entries(params ?? {}).map(([key, value]) => [
      title_case(key),
      arr_to_str(value),
    ])
</script>

<h1 style="font-size: 2.5em">{dataset.name}</h1>

<section class="meta-info">
  {#if dataset.version}
    <span>Version: {dataset.version}</span>
  {/if}

  <span title={created_ago} {@attach tooltip()}>
    <Icon icon={Calendar} /> Created: {format_date(dataset.date_created)}
  </span>

  {#if dataset.date_added}
    <span title={added_ago} {@attach tooltip()}>
      <Icon icon={CalendarPlus} /> Added: {format_date(dataset.date_added)}
    </span>
  {/if}

  <span title={dataset.n_structures.toLocaleString()} {@attach tooltip()}>
    <Icon icon={Database} />
    {format_num(dataset.n_structures, `.3~s`)} structures
  </span>

  {#if dataset.n_materials}
    <span title={dataset.n_materials.toLocaleString()} {@attach tooltip()}>
      <Icon icon={Lattice} />
      {format_num(dataset.n_materials, `.3~s`)} materials
    </span>
  {/if}

  <span
    title="The dataset is {dataset.open ? `freely ` : `in`}accessible"
    {@attach tooltip()}
  >
    <Icon icon={dataset.open ? Unlock : Lock} />
    {dataset.open ? `Open` : `Closed`}
  </span>

  <span><Icon icon={License} /> {dataset.license}</span>
</section>

<section class="links">
  {#each dataset_links as [href, label, icon, title] (label)}
    {#if href}
      <a {href} {...link_props} {title} {@attach tooltip()}>
        <Icon {icon} />
        {label}
      </a>
    {/if}
  {/each}
</section>

<section class="description">
  <h2>Description</h2>
  <p>{@html dataset.description_html}</p>
</section>

{#if dataset.temperature_range || dataset.pressure_range}
  <section class="conditions">
    <h2>Conditions</h2>
    <ul>
      {#if dataset.temperature_range}
        <li>
          Temperature Range: <strong>{dataset.temperature_range}</strong>
        </li>
      {/if}
      {#if dataset.pressure_range}
        <li>
          Pressure Range: <strong>{dataset.pressure_range}</strong>
        </li>
      {/if}
    </ul>
  </section>
{/if}

{#if dataset.contains}
  <section class="derived-from">
    <h2>Derived From</h2>
    <ol>
      {#each dataset.contains as source (source)}
        {@const contained_data = DATASETS[source]}
        <li>
          <a href="/data/{contained_data.slug}">{contained_data.name}</a>
        </li>
      {/each}
    </ol>
  </section>
{/if}

{#if dataset.method}
  <section class="method-info">
    <h2>Methodology</h2>
    <ul>
      <li>
        Method: <strong>{arr_to_str(dataset.method)}</strong>
      </li>
      {#if dataset.params}
        {#each format_params(dataset.params) as [key, value] (key)}
          <li>
            {key}: <strong>{value}</strong>
          </li>
        {/each}
      {/if}
    </ul>
  </section>
{/if}

{#if dataset.created_by && dataset.created_by.length > 0}
  <section>
    <h2>Authors</h2>
    <ol>
      {#each dataset.created_by as person (person.name)}
        <li>
          <span>{person.name}</span>
          {#if person.affiliation}
            <span class="affiliation">({person.affiliation})</span>
          {/if}
          {#if person.email}
            <a href="mailto:{person.email}" aria-label="Email">
              <Icon icon={Email} />
            </a>
          {/if}
          {#if person.github}
            <a href={person.github} {...link_props} aria-label="GitHub">
              <Icon icon={GitHub} />
            </a>
          {/if}
          {#if person.orcid}
            <a href={person.orcid} {...link_props} aria-label="ORCID">
              <Icon icon={ORCID} />
            </a>{/if}
          {#if person.url}
            <a href={person.url} {...link_props} aria-label="Website">
              <Icon icon={Globe} />
            </a>{/if}
        </li>
      {/each}
    </ol>
  </section>
{/if}

{#if dataset.slug === `mptrj`}
  <MPtrjTargetDistros />
{/if}

<p>
  See incorrect or missing data? Suggest an edit to
  <a href="{pkg.repository}/blob/-/data/datasets.yml" {...link_props}> datasets.yml </a>
</p>

<style>
  h2 {
    margin: 1em auto 0;
  }
  .meta-info,
  .links {
    display: flex;
    flex-wrap: wrap;
    gap: 3ex;
    place-content: center;
    margin: 2em auto;
  }
  section.method-info ul {
    display: flex;
    flex-wrap: wrap;
    gap: 1em;
    padding: 0;
    list-style: none;
  }
  section.method-info ul li {
    background-color: var(--nav-bg);
    padding: 2pt 6pt;
    border-radius: 3pt;
    text-align: center;
    margin: 0;
    font-weight: lighter;
    max-width: 12em;
  }
  .links a {
    padding: 0 5pt;
    background-color: var(--nav-bg);
    border-radius: 5px;
    color: var(--text-color);
  }
  .affiliation {
    color: gray;
    font-weight: lighter;
  }
  ul li {
    overflow: hidden;
    white-space: nowrap;
    text-overflow: ellipsis;
  }
</style>
