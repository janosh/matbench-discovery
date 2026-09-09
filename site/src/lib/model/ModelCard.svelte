<script lang="ts">
  import type { Label, ModelData } from '$lib/types'
  import DATASETS from '$data/datasets.yml'
  import { get_nested_number, label_data_path } from '$lib/metrics'
  import { ACTIVE_MODELS } from '$lib/models.svelte'
  import { model_metric_ranks, rank_color, RANKED_METRICS } from '$lib/rankings'
  import pkg from '$site/package.json'
  import { format_num } from 'matterviz/labels'
  import { Icon, Popover } from 'svelte-widgets'
  import {
    Calendar,
    CalendarCheck,
    Database,
    Directory,
    Docs,
    Download,
    Forest,
    GitHub,
    Info,
    NeuralNetwork,
    Paper,
  } from 'svelte-widgets/icons'
  import { tooltip } from 'svelte-widgets/attachments'

  let {
    model,
    metrics,
    sort_by,
    title_style = ``,
  }: {
    model: ModelData
    metrics: readonly Label[]
    sort_by: string // metric label key (or `Model`) highlighted in the metrics list
    title_style?: string
  } = $props()

  let { model_name, model_key, model_params, training_sets } = $derived(model)
  let ranks = $derived(model_metric_ranks(model_key, ACTIVE_MODELS, RANKED_METRICS))

  let links = $derived([
    [model.repo, `Repo`, GitHub],
    [model.paper, `Paper`, Paper],
    [model.docs, `Docs`, Docs],
    [model.checkpoint_url, `Checkpoint`, Download],
    [`${pkg.repository}/tree/HEAD/models/${model.dirname}`, `Files`, Directory],
  ] as const)
  let n_model_params = $derived(format_num(model_params, `.3~s`))
</script>

<h2 id={model_key} style={title_style}>
  <a href="/models/{model_key}">{model_name}</a>
</h2>
<nav>
  {#each links.filter( ([href]) => href?.startsWith(`http`) ) as [href, title, link_icon] (title)}
    <a {href} target="_blank" rel="noopener">
      <Icon icon={link_icon} />
      {title}
    </a>
  {/each}
</nav>

<section class="metadata">
  <span style="grid-column: span 2">
    <Icon icon={Database} />
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
    <Icon icon={Calendar} />
    Added {model.dates.benchmark_added}
  </span>
  {#if model.dates.paper_published}
    <span title="Date published">
      <Icon icon={CalendarCheck} />
      Published {model.dates.paper_published}
    </span>
  {/if}
  <span>
    <Icon icon={NeuralNetwork} />
    {n_model_params} params
  </span>
  {#if (model.n_estimators ?? 1) > 1}
    <span>
      <Icon icon={Forest} />
      Ensemble of {model.n_estimators}
      <span
        title="This result used a model ensemble with {model.n_estimators} members with {n_model_params} parameters each."
        {@attach tooltip()}
      >
        &nbsp;<Icon icon={Info} />
      </span>
    </span>
  {/if}
</section>

<section class="metrics">
  <h3 style="margin: 0; font-weight: normal">Metrics</h3>
  <ul>
    {#each metrics as metric (metric.key)}
      {@const { key, label, unit, description } = metric}
      <!-- resolve by the label's own data path so any metric works (RMSD lives under
      metrics.geo_opt.symprec=1e-2, which a hardcoded section merge would miss) -->
      {@const value = get_nested_number(model, label_data_path(metric))}
      {@const rank_entry = ranks.find((entry) => entry.metric.key === key)}
      <li class:active={sort_by == key}>
        <Popover trigger_mode="hover" trap_focus={false} aria-label="Metric description">
          {#snippet trigger(trigger_props)}
            <span
              style="display: flex; width: 100%; justify-content: space-between"
              role="button"
              tabindex="0"
              {...trigger_props}
              ><label for={key}>{@html label}</label>
              <strong>
                {#if value === undefined || isNaN(value)}
                  n/a
                {:else}
                  {format_num(value)}
                  <small>{unit ?? ``}</small>
                {/if}
              </strong></span
            >
          {/snippet}
          {@html description}
        </Popover>
        {#if value !== undefined && !isNaN(value) && rank_entry}
          <a
            class="metric-rank"
            href={rank_entry.metric.rank_href}
            aria-label="{key}: rank {rank_entry.rank} of {rank_entry.n_models}"
            ><b style:color={rank_color(rank_entry.rank, rank_entry.n_models)}
              >#{rank_entry.rank}</b
            ><span>/{rank_entry.n_models}</span></a
          >
        {/if}
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
  small {
    font-weight: 100;
    font-size: 8pt;
  }
  section.metrics > ul {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(11em, 1fr));
    justify-content: space-between;
    gap: 3pt 1em;
    list-style: none;
    padding: 0;
  }
  section.metrics > ul > li {
    font-weight: lighter;
    display: flex;
    justify-content: space-between;
    align-items: baseline;
    gap: 3pt;
  }
  section.metrics > ul > li :is(label, strong) {
    padding: 0 4pt;
    border-radius: 3pt;
  }
  section.metrics > ul > li.active label {
    font-weight: bold;
  }
  .metric-rank {
    font-size: 0.75em;
    white-space: nowrap;
    color: var(--text-secondary);
  }
  /* keep long unbroken words from widening the ModelCard container */
  :is(section, nav) {
    word-break: break-word;
  }
</style>
