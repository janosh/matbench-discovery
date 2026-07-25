<script lang="ts">
  import type { Label, ModelData } from '$lib'
  import { AuthorBrief, DATASETS } from '$lib'
  import { parse_dependency_spec } from '$lib/environment'
  import { get_nested_number, label_data_path } from '$lib/metrics'
  import pkg from '$site/package.json'
  import { format_num, Icon } from 'matterviz'
  import { tooltip } from 'svelte-multiselect/attachments'
  import type { HTMLAttributes } from 'svelte/elements'
  import { fade, slide } from 'svelte/transition'

  let {
    model,
    metrics,
    sort_by,
    show_details = $bindable(false),
    metrics_style = ``,
    title_style = ``,
    ...rest
  }: HTMLAttributes<HTMLElementTagNameMap[`section`]> & {
    model: ModelData
    metrics: readonly Label[]
    sort_by: keyof ModelData
    show_details?: boolean
    metrics_style?: string
    title_style?: string
  } = $props()

  let { model_name, model_key, model_params, training_sets } = $derived(model)
  let env_packages = $derived(model.environment.dependencies.map(parse_dependency_spec))

  let links = $derived([
    [model.repo, `Repo`, `GitHub`],
    [model.paper, `Paper`, `Paper`],
    [model.docs, `Docs`, `Docs`],
    [model.checkpoint_url, `Checkpoint`, `Download`],
    [`${pkg.repository}/tree/HEAD/models/${model.dirname}`, `Files`, `Directory`],
  ] as const)
  const target = { target: `_blank`, rel: `noopener` }
  let n_model_params = $derived(format_num(model_params, `.3~s`))
  let expand_title = $derived(
    `${show_details ? `Hide` : `Show`} authors and package versions`,
  )
</script>

<h2 style={title_style}>
  <a href="/models/{model_key}">{model_name}</a>
  <button
    onclick={() => (show_details = !show_details)}
    aria-expanded={show_details}
    aria-label={expand_title}
    title={expand_title}
    style={title_style}
  >
    <Icon icon="Arrow{show_details ? `Up` : `Down`}" />
  </button>
</h2>
<nav>
  {#each links.filter( ([href]) => href?.startsWith(`http`) ) as [href, title, icon] (title)}
    <a {href} {...target}>
      <Icon {icon} />
      {title}
    </a>
  {/each}
</nav>

<section class="metadata" {...rest}>
  <span style="grid-column: span 2">
    <Icon icon="Database" />
    Training data:
    {#each training_sets as train_set_key, idx (train_set_key)}
      {#if idx > 0}
        &nbsp;+&nbsp;
      {/if}
      {@const { n_structures, name, slug, n_materials } = DATASETS[train_set_key]}
      {@const n_mat_str = n_materials ? ` from ${format_num(n_materials)} materials` : ``}
      <a
        href="/data/{slug}"
        title="{name}: {format_num(n_structures)} structures{n_mat_str}"
        {@attach tooltip()}
      >
        {train_set_key}
      </a>
    {/each}
  </span>
  <span title="Date added">
    <Icon icon="Calendar" />
    Added {model.dates.benchmark_added}
  </span>
  {#if model.dates.paper_published}
    <span title="Date published">
      <Icon icon="CalendarCheck" />
      Published {model.dates.paper_published}
    </span>
  {/if}
  <span>
    <Icon icon="NeuralNetwork" />
    {n_model_params} params
  </span>
  {#if (model.n_estimators ?? 1) > 1}
    <span>
      <Icon icon="Forest" />
      Ensemble of {model.n_estimators}
      <span
        title="This result used a model ensemble with {model.n_estimators} members with {n_model_params} parameters each."
        {@attach tooltip()}
      >
        &nbsp;<Icon icon="Info" />
      </span>
    </span>
  {/if}
</section>
{#if show_details}
  <div transition:fade={{ duration: 200 }}>
    <section transition:slide={{ duration: 200 }}>
      <h3>Authors</h3>
      <ul>
        {#each model.authors as author (author.name)}
          <li><AuthorBrief {author} /></li>
        {/each}
      </ul>
      {#if model.trained_by}
        <h3>Trained By</h3>
        <ul>
          {#each model.trained_by as author (author.name)}
            <li><AuthorBrief {author} /></li>
          {/each}
        </ul>
      {/if}
    </section>
    <section transition:slide={{ duration: 200 }}>
      <h3>Package versions</h3>
      <ul>
        {#each env_packages as { name, detail, href } (name + detail)}
          <li style="font-size: smaller">
            {name}{#if detail}:
              <a {href} {...target}>{detail}</a>{/if}
          </li>
        {/each}
      </ul>
    </section>
  </div>
{/if}

<section class="metrics" style={metrics_style || null}>
  <h3 style="margin: 0; font-weight: normal">Metrics</h3>
  <ul>
    {#each metrics as metric (metric.key)}
      {@const { key, label, unit, description } = metric}
      <!-- resolve by the label's own data path so any metric works (RMSD lives under
      metrics.geo_opt.symprec=1e-2, which a hardcoded section merge would miss) -->
      {@const value = get_nested_number(model, label_data_path(metric))}
      <li
        class:active={sort_by == key}
        title={description}
        {@attach tooltip({ allow_html: true })}
      >
        <label for={key}>{@html label ?? key}</label>
        <strong>
          {#if value === undefined || isNaN(value)}
            n/a
          {:else}
            {format_num(value)}
            <small>{unit ?? ``}</small>
          {/if}
        </strong>
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
  nav > a {
    display: inline-flex;
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
    display: grid;
    /* cap track width so label and value stay close (widest current pair is ~9.2em:
    'MAE 0.044 eV / atom'); leftover space goes between columns (justify-content),
    not inside cells where it'd stretch the label-value gap */
    grid-template-columns: repeat(auto-fill, minmax(9em, 9.5em));
    justify-content: space-between;
    gap: 3pt 1em;
    list-style: none;
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
  section.metrics > ul > li.active > label {
    font-weight: bold;
  }
  /* prevent long from increasing ModelCard container width */
  :is(section, div, nav) {
    word-break: break-word;
  }
</style>
