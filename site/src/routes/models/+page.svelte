<script lang="ts">
  import { repository } from '$site/package.json'
  import Icon from '@iconify/svelte'
  import type { PageData } from './$types'

  export let data: PageData
</script>

<h1>Models</h1>

<ol>
  {#each data.model_metas as [key, meta], idx}
    {@const { repo, doi, preprint, url, date_added } = meta}
    <li>
      <h2>{meta.model_name}</h2>
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
        &nbsp;&bull;&nbsp; Benchmark version: {meta.matbench_discovery_version}
      </p>
      <strong>Authors</strong>
      <section>
        <ul>
          {#each meta.authors as { name, email, orcid, affiliation, url }}
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
        <strong>Package versions</strong>
        <ul>
          {#each Object.entries(meta.requirements) as [name, version]}
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
      <!-- TODO add table with performance metrics (F1, Acc, Recall, Precision) for each model -->
    </li>
  {/each}
</ol>

<style>
  ol {
    display: grid;
    gap: 2em;
    list-style: none;
  }
  ol > li {
    background-color: rgba(255, 255, 255, 0.05);
    padding: 3pt 9pt 5pt;
  }
  ol > li > h2 {
    margin: 5pt 0 1ex;
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
    gap: 0.5em;
    place-items: center;
  }
  strong {
    display: block;
    margin: 1em 0 5pt;
  }
</style>
