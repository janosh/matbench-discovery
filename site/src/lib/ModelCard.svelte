<script lang="ts">
  import { repository } from '$site/package.json'
  import Icon from '@iconify/svelte'
  import { fade, slide } from 'svelte/transition'
  import { pretty_num } from 'sveriodic-table/labels'
  import type { ModelData, ModelStatLabel } from '.'

  export let key: string
  export let data: ModelData
  export let stats: ModelStatLabel[] // [key, label, unit][]
  export let sort_by: keyof ModelData
  export let show_details: boolean = false

  $: ({ model_name, repo, doi, preprint, url, date_added } = data)
  $: ({ missing_preds, missing_percent, hyperparams, notes } = data)

  const links = [
    [repo, `Repo`, `octicon:mark-github`],
    [doi, `DOI`, `academicons:doi`],
    [preprint, `Preprint`, `ion:ios-paper`],
    [url, `Website`, `ion:ios-globe`],
    [`${repository}/blob/-/models/${key}`, `Files`, `octicon:file-directory`],
  ]
  const target = { target: `_blank`, rel: `noopener` }
</script>

<h2>
  {model_name}
  <button
    on:click={() => (show_details = !show_details)}
    title="{show_details ? 'Hide' : 'Show'} authors and package versions"
  >
    <!-- change between expand/collapse icon -->
    <Icon icon={show_details ? `ion:ios-arrow-up` : `ion:ios-arrow-down`} inline />
  </button>
</h2>
<nav>
  {#each links as [href, title, icon]}
    <span>
      <Icon {icon} inline />
      <a {href} {...target}>{title}</a>
    </span>
  {/each}
</nav>
<p>
  Date added: {date_added}
  &nbsp;&bull;&nbsp; Benchmark version: {data.matbench_discovery_version}
  &nbsp;&bull;&nbsp; Missing predictions:
  {pretty_num(missing_preds)}
  <small>({missing_percent})</small>
</p>
{#if show_details}
  <div transition:fade|fly={{ duration: 200 }}>
    <section transition:slide={{ duration: 200 }}>
      <h3 class="toc-exclude">Authors</h3>
      <ul>
        {#each data.authors as { name, email, orcid, affiliation, url }}
          <li>
            <span title={affiliation}>{name}</span>
            {#if email}
              [<a href="mailto:{email}">email</a>]
            {/if}
            {#if orcid}
              [<a href={orcid}>Orcid</a>]
            {/if}
            {#if url}
              [<a href={url}>web</a>]
            {/if}
          </li>
        {/each}
      </ul>
    </section>
    <section transition:slide={{ duration: 200 }}>
      <h3 class="toc-exclude">Package versions</h3>
      <ul>
        {#each Object.entries(data.requirements ?? {}) as [name, version]}
          <li>
            {#if ![`aviary`].includes(name)}
              {@const href = `https://pypi.org/project/${name}/${version}`}
              {name}: <a {href} {...target}>{version}</a>
            {:else}
              {name}: {version}
            {/if}
          </li>
        {/each}
      </ul>
    </section>
  </div>
{/if}
<section class="metrics">
  <h3 class="toc-exclude">Metrics</h3>
  <ul>
    {#each stats as { key, label, unit }}
      <li class:active={sort_by == key}>
        {@html label ?? key} = {data[key]}
        {unit ?? ``}
      </li>
    {/each}
  </ul>
</section>
{#if hyperparams && show_details}
  <section>
    <h3 class="toc-exclude">Hyperparameters</h3>
    <ul>
      {#each Object.entries(hyperparams) as [key, value]}
        <li>
          {key}:
          {#if typeof value == `object`}
            <ul>
              {#each Object.entries(value) as [k, v]}
                <li><code>{k} = {v}</code></li>
              {/each}
            </ul>
          {:else}
            {value}
          {/if}
        </li>
      {/each}
    </ul>
  </section>
{/if}
{#if notes && show_details}
  <section>
    <h3 class="toc-exclude">Notes</h3>
    <ul>
      {#each Object.values(notes) as note}
        <li>{@html note}</li>
      {/each}
    </ul>
  </section>
{/if}

<style>
  h2 {
    margin: 8pt 0 1em;
    text-align: center;
  }
  button {
    background: none;
    padding: 0;
    font: inherit;
  }
  h3 {
    margin: 1em 0 0;
  }
  div h3 {
    margin: 0;
  }
  ul {
    list-style: disc;
  }
  nav {
    font-weight: 250;
    display: flex;
    gap: 5pt 1em;
    flex-wrap: wrap;
    place-content: center;
  }
  nav > span {
    display: flex;
    gap: 6pt;
    place-items: center;
  }
  div {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 15pt;
    justify-content: space-between;
  }
  small {
    font-weight: lighter;
  }
  section.metrics > ul {
    display: flex;
    flex-wrap: wrap;
    gap: 0 1em;
    list-style: none;
  }
  section.metrics > ul > li.active {
    font-weight: 700;
  }
</style>
