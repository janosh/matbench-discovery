<script lang="ts">
  import type { ModelData, ModelStatLabel } from '$lib'
  import { AuthorBrief, DATASETS } from '$lib'
  import { repository } from '$site/package.json'
  import { pretty_num } from 'elementari'
  import 'iconify-icon'
  import { Tooltip } from 'svelte-zoo'
  import { fade, slide } from 'svelte/transition'

  interface Props {
    model: ModelData
    stats: ModelStatLabel[]
    sort_by: keyof ModelData
    show_details?: boolean
    style?: string | null
    metrics_style?: string | null
  }
  let {
    model,
    stats,
    sort_by,
    show_details = $bindable(false),
    style = null,
    metrics_style = null,
  }: Props = $props()

  let { model_name, model_key, model_params } = $derived(model)
  let { notes = {}, training_set, n_estimators } = $derived(model)
  let all_metrics = $derived({
    ...(model.metrics?.discovery?.full_test_set ?? {}),
    ...(typeof model.metrics?.phonons == `object`
      ? model.metrics?.phonons.kappa_103
      : {}),
  })
  let { missing_preds, missing_percent } = $derived(all_metrics)

  let links = $derived([
    [model.repo, `Repo`, `octicon:mark-github`],
    [model.paper, `Paper`, `ion:ios-paper`],
    [model.url, `Docs`, `ion:ios-globe`],
    [model.checkpoint_url, `Checkpoint`, `mdi:download-box`],
    [`${repository}/blob/-/models/${model.dirname}`, `Files`, `octicon:file-directory`],
  ])
  const target = { target: `_blank`, rel: `noopener` }
  let n_model_params = $derived(pretty_num(model_params, `.3~s`))
</script>

<h2 id={model_key} {style}>
  <a href="/models/{model_key}">{model_name}</a>
  <button
    onclick={() => (show_details = !show_details)}
    title="{show_details ? `Hide` : `Show`} authors and package versions"
  >
    <!-- change between expand/collapse icon -->
    <iconify-icon icon={show_details ? `ion:ios-arrow-up` : `ion:ios-arrow-down`} inline
    ></iconify-icon>
  </button>
</h2>
<nav>
  {#each links.filter( ([href]) => href?.startsWith(`http`), ) as [href, title, icon] (title)}
    <span>
      <iconify-icon {icon} inline></iconify-icon>
      <a {href} {...target}>{title}</a>
    </span>
  {/each}
</nav>

<section class="metadata">
  {#if training_set}
    <span style="grid-column: span 2;">
      <iconify-icon icon="mdi:database" inline></iconify-icon>
      Training data:
      {#each training_set as train_set_key, idx (train_set_key)}
        {#if idx > 0}
          &nbsp;+&nbsp;
        {/if}
        {@const training_set = DATASETS[train_set_key]}
        {@const { n_structures, url, title, n_materials } = training_set}
        {@const pretty_n_mat =
          typeof n_materials == `number` ? pretty_num(n_materials) : n_materials}
        {@const n_mat_str = n_materials ? ` from ${pretty_n_mat} materials` : ``}
        <Tooltip text="{title}: {pretty_num(n_structures)} structures{n_mat_str}">
          <a href={url}>{train_set_key}</a>
        </Tooltip>
      {/each}
    </span>
  {/if}
  <span title="Date added">
    <iconify-icon icon="ion:ios-calendar" inline></iconify-icon>
    Added {model.date_added}
  </span>
  {#if model.date_published}
    <span title="Date published">
      <iconify-icon icon="ri:calendar-check-line" inline></iconify-icon>
      Published {model.date_published}
    </span>
  {/if}
  <span>
    <iconify-icon icon="eos-icons:neural-network" inline></iconify-icon>
    {n_model_params} params
  </span>
  {#if n_estimators > 1}
    <span>
      <iconify-icon icon="material-symbols:forest" inline></iconify-icon>
      Ensemble of {n_estimators > 1 ? `${n_estimators}` : ``}
      <Tooltip
        text="This result used a model ensemble with {n_estimators} members with {n_model_params} parameters each."
      >
        &nbsp;<iconify-icon icon="octicon:info-24" inline></iconify-icon>
      </Tooltip>
    </span>
  {/if}
  <span>
    <iconify-icon icon="fluent:missing-metadata-24-regular" inline></iconify-icon>
    {#if missing_preds}
      Missing preds:
      {pretty_num(missing_preds, `,.0f`)}
      <small>({missing_percent})</small>
    {/if}
    {#if notes.html?.missing_preds}
      <Tooltip text={notes.html.missing_preds ?? ``} tip_style="font-size: 9pt;">
        &nbsp;<iconify-icon icon="octicon:info-24" inline></iconify-icon>
      </Tooltip>
    {/if}
  </span>
</section>
{#if show_details}
  <div transition:fade|fly={{ duration: 200 }}>
    <section transition:slide={{ duration: 200 }}>
      <h3>Authors</h3>
      <ul>
        {#each model.authors as author (author.name)}
          <li>
            <AuthorBrief {author} />
          </li>
        {/each}
      </ul>
      {#if model.trained_by}
        <h3>Trained By</h3>
        <ul>
          {#each model.trained_by as author (author.name)}
            <li>
              <AuthorBrief {author} />
            </li>
          {/each}
        </ul>
      {/if}
    </section>
    <section transition:slide={{ duration: 200 }}>
      <h3>Package versions</h3>
      <ul>
        {#each Object.entries(model.requirements ?? {}) as [name, version] (name)}
          {@const [href, link_text] = version.startsWith(`http`)
            ? // version.split(`/`).at(-1) assumes final part after / of URL is the package version, as is the case for GitHub releases
              [version, version.split(`/`).at(-1)]
            : [`https://pypi.org/project/${name}/${version}`, version]}
          <li style="font-size: smaller;">
            {name}: <a {href} {...target}>{link_text}</a>
          </li>
        {/each}
      </ul>
    </section>
  </div>
{/if}
<section class="metrics" style={metrics_style}>
  <h3>Metrics</h3>
  <ul>
    <!-- hide run time if value is 0 (i.e. not available) -->
    {#each stats as { key, label, unit } (key)}
      <li class:active={sort_by == key}>
        <label for={key}>{@html label ?? key}</label>
        <strong>{pretty_num(all_metrics[key])} <small>{unit ?? ``}</small></strong>
      </li>
    {/each}
  </ul>
</section>

<style>
  h2 {
    margin: 8pt 0 0;
    text-align: center;
    border-radius: 5pt;
  }
  h2 a {
    color: inherit;
  }
  button {
    background: none;
    padding: 0;
    font: inherit;
  }
  h3 {
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
  section.metadata {
    display: grid;
    gap: 9pt 5pt;
    grid-template-columns: 1fr 1fr;
    font-size: 0.95em;
    align-content: center;
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
    padding: 0;
  }
  section.metrics > ul > li {
    font-weight: lighter;
    display: flex;
    justify-content: space-between;
  }
  section.metrics > ul > li > :is(label, strong) {
    padding: 0 4pt;
    border-radius: 3pt;
  }
  section.metrics > ul > li > strong {
    background-color: rgba(0, 0, 0, 0.25);
  }
  section.metrics > ul > li.active > strong {
    background-color: teal;
  }
  section.metrics > ul > li.active > label {
    font-weight: bold;
  }
  /* prevent long from increasing ModelCard container width */
  :is(section, div, nav) {
    word-break: break-word;
  }
</style>
