<script lang="ts">
  import DATASETS from '$data/datasets.yml'
  import { arr_to_str, format_date } from '$lib'
  import { format_relative_time, title_case } from '$lib/labels'
  import { format_num } from 'matterviz/labels'
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
    Key,
    Lattice,
    License,
    ORCID,
    Tag,
  } from 'svelte-widgets/icons'
  import pkg from '$site/package.json'
  import { tooltip } from 'svelte-widgets/attachments'
  import type { PageData } from './$types'
  import MPtrjTargetDistros from './MPtrjTargetDistros.svelte'
  import WbmDetails from './WbmDetails.svelte'

  let { data }: { data: PageData } = $props()
  let dataset = $derived(data.dataset)
  const link_props = { target: `_blank`, rel: `noopener noreferrer` }
  const source_url = `${pkg.repository}/blob/main/data/datasets.yml`

  const metadata = $derived([
    [Tag, dataset.version ? `Version: ${dataset.version}` : null],
    [
      Calendar,
      `Created: ${format_date(dataset.date_created)}`,
      format_relative_time(dataset.date_created),
    ],
    [
      CalendarPlus,
      dataset.date_added ? `Added: ${format_date(dataset.date_added)}` : null,
      format_relative_time(dataset.date_added),
    ],
    [
      Database,
      dataset.n_structures === null
        ? `Unknown number of structures`
        : `${format_num(dataset.n_structures, `.3~s`)} structures`,
      dataset.n_structures?.toLocaleString() ?? `Structure count not reported`,
    ],
    [
      Lattice,
      dataset.n_materials ? `${format_num(dataset.n_materials, `.3~s`)} materials` : null,
      dataset.n_materials?.toLocaleString(),
    ],
    [Key, `${title_case(dataset.access)} access`],
    [License, dataset.license],
  ] as const)
  let dataset_links = $derived([
    [dataset.url, `Website`, Globe, `View dataset website`],
    [dataset.download_url, `Download`, Download, `Download dataset`],
    [dataset.doi, `DOI`, DOI, `Digital Object Identifier`],
    [source_url, `Source`, Code, `View source YAML file`],
  ] as const)
</script>

<h1 style="font-size: 2.5em">{dataset.name}</h1>

<section class="meta-info">
  {#each metadata as [icon, label, title] (icon)}
    {#if label !== null}
      <span {title} {@attach title ? tooltip() : undefined}>
        <Icon {icon} />
        {label}
      </span>
    {/if}
  {/each}
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
  <h2 id="description">Description</h2>
  {@html dataset.description_html}
</section>

{#each Object.entries(dataset.notes_html ?? {}) as [title, note] (title)}
  <details>
    <summary>{title}</summary>
    {@html note}
  </details>
{/each}

{#if dataset.temperature_range || dataset.pressure_range}
  <section class="conditions">
    <h2 id="conditions">Conditions</h2>
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
    <h2 id="derived-from">Derived From</h2>
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
    <h2 id="methodology">Methodology</h2>
    <ul>
      <li>
        Method: <strong>{arr_to_str(dataset.method)}</strong>
      </li>
      {#each Object.entries(dataset.params ?? {}) as [key, value] (key)}
        <li>{title_case(key)}: <strong>{arr_to_str(value)}</strong></li>
      {/each}
    </ul>
  </section>
{/if}

{#if dataset.created_by && dataset.created_by.length > 0}
  <section>
    <h2 id="authors">Authors</h2>
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
{:else if dataset.slug === `wbm`}
  <WbmDetails />
{/if}

<p>
  See incorrect or missing data? Suggest an edit to
  <a href={source_url} {...link_props}>datasets.yml</a>
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
    :global(svg) {
      width: 1.2em;
      transform: translateY(-2px);
    }
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
