<script lang="ts">
  import { repository } from '$site/package.json'
  import Icon from '@iconify/svelte'
  import { pretty_num } from 'sveriodic-table/labels'
  import type { ModelMetadata } from './types'

  export let key: string
  export let data: ModelMetadata

  const { model_name, repo, doi, preprint, url, date_added } = data
  const { missing_preds, test_set_size } = data
</script>

<h2>{model_name}</h2>
<nav>
  {#each [[repo, `Repo`, `octicon:mark-github`], [`${repository}/tree/main/models/${key}`, `Submission files`, `octicon:file-directory`], [doi, `DOI`, `academicons:doi`], [preprint, `Preprint`, `ion:ios-paper`], [url, `Website`, `ion:ios-globe`]] as [href, title, icon]}
    <span>
      <Icon {icon} inline />
      <a {href}>{title}</a>
    </span>
  {/each}
</nav>
<p>
  Date added: {new Date(date_added).toISOString().split(`T`)[0]}
  &nbsp;&bull;&nbsp; Benchmark version: {data.matbench_discovery_version}
  &nbsp;&bull;&nbsp; Missing predictions:
  {pretty_num(missing_preds)}
  <small>({((100 * missing_preds) / test_set_size).toFixed(2)}%)</small>
</p>
<div>
  <section>
    <h3>Authors</h3>
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
    <h3>Package versions</h3>
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

<!-- TODO add table with performance metrics (F1, Acc, Recall, Precision) for each model -->
<style>
  h2 {
    margin: 5pt 0 1ex;
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
    display: flex;
    gap: 15pt;
    margin: 1em 0;
    justify-content: space-between;
  }
  small {
    font-weight: lighter;
  }
</style>
