<script lang="ts">
  import { arr_to_str, calculate_days_ago, format_date, Icon, slugify } from '$lib'
  import type { Dataset } from '$lib/types'
  import pkg from '$site/package.json'
  import { format_num } from 'matterviz'
  import { Tooltip } from 'svelte-zoo'
  import { titles_as_tooltips } from 'svelte-zoo/actions'
  import MptrjTargetDistros from './MptrjTargetDistros.svelte'

  interface Props {
    data: { dataset: Dataset }
  }
  let { data }: Props = $props()
  const { dataset } = data
  const ext_link_props = { target: `_blank`, rel: `noopener noreferrer` }

  let days_created = calculate_days_ago(dataset.date_created)
  let days_added = calculate_days_ago(dataset.date_added ?? ``)

  // Format the params object into a readable list
  function format_params(params: Record<string, unknown> | undefined): string[] {
    if (!params) return []

    return Object.entries(params).map(([key, value]) => {
      const formatted_key = key
        .split(`_`)
        .map((word) => word.charAt(0).toUpperCase() + word.slice(1))
        .join(` `)
      return `${formatted_key}: ${arr_to_str(value)}`
    })
  }
</script>

<h1 style="font-size: 2.5em">{dataset.name}</h1>

<section class="meta-info">
  {#if dataset.version}
    <div>Version: {dataset.version}</div>
  {/if}

  <div>
    <Icon icon="Calendar" />
    Created: <Tooltip text="{days_created} days ago">
      {format_date(dataset.date_created)}
    </Tooltip>
  </div>

  {#if dataset.date_added}
    <div>
      <Icon icon="CalendarPlus" />
      Added: <Tooltip text="{days_added} days ago">
        {format_date(dataset.date_added)}
      </Tooltip>
    </div>
  {/if}

  <div>
    <Icon icon="Database" />
    <Tooltip text={dataset.n_structures.toLocaleString()}>
      {format_num(dataset.n_structures, `.3~s`)}
    </Tooltip> structures
  </div>

  {#if dataset.n_materials}
    <div>
      <Icon icon="Lattice" />
      <Tooltip text={dataset.n_materials.toLocaleString()}>
        {format_num(dataset.n_materials, `.3~s`)}
      </Tooltip> materials
    </div>
  {/if}

  <div>
    <Icon icon={dataset.open ? `Unlock` : `Lock`} />
    {dataset.open ? `Open` : `Closed`}
  </div>

  <div>
    <Icon icon="License" />
    {dataset.license}
  </div>
</section>

<section class="links">
  <a
    href={dataset.url}
    {...ext_link_props}
    title="View dataset website"
    use:titles_as_tooltips
  >
    <Icon icon="Globe" /> Website
  </a>

  {#if dataset.download_url}
    <a
      href={dataset.download_url}
      {...ext_link_props}
      title="Download dataset"
      use:titles_as_tooltips
    >
      <Icon icon="Download" /> Download
    </a>
  {/if}

  {#if dataset.doi}
    <a
      href={dataset.doi}
      {...ext_link_props}
      title="Digital Object Identifier"
      use:titles_as_tooltips
    >
      <Icon icon="DOI" /> DOI
    </a>
  {/if}

  <a
    href="{pkg.repository}/blob/main/data/datasets.yml"
    {...ext_link_props}
    title="View source YAML file"
    use:titles_as_tooltips
  >
    <Icon icon="Code" /> Source
  </a>
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

{#if dataset.derived_from}
  <section class="derived-from">
    <h2>Derived From</h2>
    <ul>
      {#each dataset.derived_from as source (source)}
        <li>
          <a href="/data/{slugify(source)}">{source}</a>
        </li>
      {/each}
    </ul>
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
        {#each format_params(dataset.params) as param (param)}
          {@const [key, value] = param.split(`:`)}
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
              <Icon icon="Contact" />
            </a>
          {/if}
          {#if person.github}
            <a href={person.github} {...ext_link_props} aria-label="GitHub">
              <Icon icon="GitHub" />
            </a>
          {/if}
          {#if person.orcid}
            <a href={person.orcid} {...ext_link_props} aria-label="ORCID">
              <Icon icon="Orcid" />
            </a>{/if}
          {#if person.url}
            <a href={person.url} {...ext_link_props} aria-label="Website">
              <Icon icon="Globe" />
            </a>{/if}
        </li>
      {/each}
    </ol>
  </section>
{/if}

{#if dataset.slug === `mptrj`}
  <MptrjTargetDistros />
{/if}

<p>
  See incorrect or missing data? Suggest an edit to
  <a href="{pkg.repository}/blob/-/data/datasets.yml" {...ext_link_props}>
    datasets.yml
  </a>
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
  section:is(.method-info) ul {
    display: flex;
    flex-wrap: wrap;
    gap: 1em;
    padding: 0;
    list-style: none;
  }
  section:is(.method-info) ul li {
    background-color: rgba(255, 255, 255, 0.1);
    padding: 2pt 6pt;
    border-radius: 3pt;
    text-align: center;
    margin: 0;
    font-weight: lighter;
    max-width: 12em;
  }
  .links a {
    gap: 5px;
    padding: 2pt 6pt;
    background-color: rgba(255, 255, 255, 0.1);
    border-radius: 5px;
    color: lightgray;
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
