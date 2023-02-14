<script lang="ts">
  import { repository } from '$site/package.json'
  import Icon from '@iconify/svelte'
  import { pretty_num } from 'elementari/labels'
  import { fade, slide } from 'svelte/transition'
  import type { ModelData, ModelStatLabel } from '.'

  export let data: ModelData
  export let idx: number
  export let stats: ModelStatLabel[] // [key, label, unit][]
  export let sort_by: keyof ModelData
  export let show_details: boolean = false

  $: ({ model_name, missing_preds, missing_percent, hyperparams, notes } = data)

  $: links = [
    [data.repo, `Repo`, `octicon:mark-github`],
    [data.doi, `DOI`, `academicons:doi`],
    [data.preprint, `Preprint`, `ion:ios-paper`],
    [data.url, `Website`, `ion:ios-globe`],
    [`${repository}/blob/-/models/${data.dir}`, `Files`, `octicon:file-directory`],
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
<em>{idx + 1}</em>
<nav>
  {#each links.filter(([href]) => href) as [href, title, icon]}
    <span>
      <Icon {icon} inline />
      <a {href} {...target}>{title}</a>
    </span>
  {/each}
</nav>
<p>
  <span
    ><Icon icon="ion:ios-calendar" inline />
    {data.date_added}
  </span>
  <span
    ><Icon icon="carbon:version" inline />
    Benchmark version: {data.matbench_discovery_version}
  </span>
  <span>
    <Icon icon="fluent:missing-metadata-24-regular" inline />
    Missing preds:
    {pretty_num(missing_preds ?? null)}
    <small>({missing_percent})</small>
  </span>
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
              <a href="mailto:{email}"><Icon icon="ion:ios-mail" inline /></a>
            {/if}
            {#if orcid}
              <a href={orcid}><Icon icon="fa-brands:orcid" inline /></a>
            {/if}
            {#if url}
              <a href={url}><Icon icon="ion:ios-globe" inline /></a>
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
        <label for={key}>{@html label ?? key}</label>
        <strong>{data[key]} <small>{unit ?? ``}</small></strong>
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
      {#each [`description`, `training`].filter((k) => k in (notes ?? {})) as key}
        <li>{@html notes[key]}</li>
      {/each}
    </ul>
  </section>
{/if}

<style>
  h2 {
    margin: 8pt 0 1em;
    text-align: center;
  }
  em {
    font-weight: lighter;
    font-style: normal;
    position: absolute;
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
  p {
    display: flex;
    gap: 3pt 12pt;
    flex-wrap: wrap;
  }
  div {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 15pt;
    justify-content: space-between;
  }
  small {
    font-weight: 100;
    font-size: 8pt;
  }
  section.metrics > ul {
    display: flex;
    flex-wrap: wrap;
    gap: 3pt 1em;
    list-style: none;
    flex-direction: column;
    max-height: 10em;
  }
  section.metrics > h3 {
    margin: 0;
  }
  section.metrics > ul > li {
    font-weight: lighter;
    display: flex;
    justify-content: space-between;
  }
  section.metrics > ul > li > strong {
    background-color: rgba(0, 0, 0, 0.3);
    padding: 0 4pt;
    border-radius: 3pt;
  }
  section.metrics > ul > li.active > strong {
    background-color: teal;
  }
</style>
