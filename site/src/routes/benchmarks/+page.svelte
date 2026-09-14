<script lang="ts">
  import MODELING_TASKS from '$pkg/modeling-tasks.yml'
  import { benchmarks } from '$lib/benchmark/data'
  import ElementHeatmap from '$lib/benchmark/ElementHeatmap.svelte'
  import { Icon } from 'svelte-widgets'
  import DataFilesDirectDownload from '../data/data-files-direct-download.md'
</script>

<h1>Benchmarks</h1>
<p class="intro">
  Explore what each benchmark measures, the data it uses, and how models compare.
</p>
<div class="benchmark-grid bleed-1400">
  {#each Object.entries(benchmarks) as [task, { key, icon, dataset }] (task)}
    <article>
      <h2>
        <a href="/benchmarks/{task}"><Icon {icon} /> {MODELING_TASKS[key].label}</a>
      </h2>
      <figure aria-label="{dataset.name} element occurrences">
        <ElementHeatmap {dataset} />
      </figure>
      <p>{MODELING_TASKS[key].description}</p>
      <p class="dataset"><strong>{dataset.name}</strong><br />{dataset.coverage}</p>
      <a href="/benchmarks/{task}#test-set">Test set and downloads →</a>
    </article>
  {/each}
</div>
<p>
  Browse the <a href="/data/sets">dataset catalog</a> for training datasets and repositories,
  access terms, licenses, and model usage. Keep benchmark test data separate from model training
  and tuning.
</p>
<DataFilesDirectDownload />

<style>
  h1,
  .intro {
    text-align: center;
  }
  .benchmark-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(min(100%, 23rem), 1fr));
    gap: 1em;
    margin-block: 2em;
    article {
      --heatmap-width: min(60%, 13.2rem);
      display: flow-root;
      padding: 1.25em;
      border: 1px solid var(--border);
      border-radius: 8px;
      background: var(--card-bg);
    }
    h2 {
      float: left;
      width: calc(100% - var(--heatmap-width) - 1rem);
      min-width: min-content;
      max-width: 100%;
      font-size: 1.25em;
      line-height: 1.25;
    }
    h2,
    p {
      margin: 0 0 1rem;
    }
    p {
      hyphens: auto;
      overflow-wrap: anywhere;
    }
    figure + p {
      clear: left;
    }
    h2 a {
      display: inline-flex;
      align-items: center;
      gap: 0.4em;
    }
    .dataset {
      color: var(--text-secondary);
      font-size: 0.9em;
    }
    figure {
      float: right;
      width: var(--heatmap-width);
      margin: 0 0 1rem 1rem;
      --elem-tile-border-radius: 1px;
    }
  }
</style>
