<script lang="ts">
  import { repository } from '$site/package.json'
  import Icon from '@iconify/svelte'
  import { pretty_num } from 'sveriodic-table/labels'
  import type { ModelData, ModelStat } from '.'

  export let key: string
  export let data: ModelData
  export let stats: [ModelStat, string?][]
  export let sort_by: keyof ModelData

  $: ({ model_name, repo, doi, preprint, url, date_added } = data)
  $: ({ missing_preds, missing_percent } = data)

  const links = [
    [repo, `Repo`, `octicon:mark-github`],
    [doi, `DOI`, `academicons:doi`],
    [preprint, `Preprint`, `ion:ios-paper`],
    [url, `Website`, `ion:ios-globe`],
    [`${repository}/blob/-/models/${key}`, `Submission files`, `octicon:file-directory`],
  ]
  const target = { target: `_blank`, rel: `noopener` }
</script>

<h2>{model_name}</h2>
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
<div>
  <section>
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
  <section>
    <h3 class="toc-exclude">Package versions</h3>
    <ul>
      {#each Object.entries(data.requirements) as [name, version]}
        <li>
          {#if ![`aviary`].includes(name)}
            {@const href = `https://pypi.org/project/${name}/${version}`}
            {name}: <a {href}>{version}</a>
          {:else}
            {name}: {version}
          {/if}
        </li>
      {/each}
    </ul>
  </section>
</div>
<section class="metrics">
  <h3 class="toc-exclude">Metrics</h3>
  <ul>
    {#each stats as [key, label]}
      <li class:active={sort_by == key}>
        {@html label ?? key} = {data[key]}
      </li>
    {/each}
  </ul>
</section>

<!-- TODO add table with performance metrics (F1, Acc, Recall, Precision) for each model -->
<style>
  h2 {
    margin: 8pt 0 1em;
    text-align: center;
  }
  h3 {
    margin: 0;
  }
  ul {
    list-style-type: disc;
  }
  nav {
    font-weight: 250;
    display: flex;
    gap: 5pt 1em;
    flex-wrap: wrap;
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
    margin: 1em 0;
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
    font-weight: 600;
    color: lightseagreen;
  }
</style>
