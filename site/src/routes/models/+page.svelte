<script lang="ts">
  import RunTimePie from '$figs/2023-01-26-model-run-times-pie.svelte'
  import type { ModelStatLabel, ModelStats } from '$lib'
  import { ModelCard } from '$lib'
  import { RadioButtons } from 'svelte-zoo'
  import { flip } from 'svelte/animate'
  import { fade } from 'svelte/transition'
  import type { PageData } from './$types'

  export let data: PageData

  let sort_by: keyof ModelStats | 'model_name' = `model_name`
  let selected = `asc`
  $: sort_factor = selected == `asc` ? -1 : 1

  $: models = data.models.sort(([_k1, m1], [_k2, m2]) => {
    if (typeof m1[sort_by] == `string`) {
      return sort_factor * -m1[sort_by].localeCompare(m2[sort_by])
    } else if (typeof m1[sort_by] == `number`) {
      return sort_factor * (m2[sort_by] - m1[sort_by])
    } else {
      console.error(`Sorting by key ${sort_by} gives unknown type: ${typeof m1[sort_by]}`)
    }
  })
  const stats: ModelStatLabel[] = [
    // key, label, unit
    [`MAE`, null, `eV / atom`],
    [`RMSE`, null, `eV / atom`],
    [`R2`, `R<sup>2</sup>`],
    [`Precision`],
    [`Recall`],
    [`F1`],
    [`date_added`, `Date added`],
    [`run_time_h`, `Run time`, `h`],
  ]
</script>

<div class="pull-left">
  <h1>Models</h1>

  <span>
    Sort <RadioButtons bind:selected options={[`asc`, `desc`]} /> by:
  </span>
  <ul>
    {#each [[`model_name`, `Model Name`], ...stats] as [key, label]}
      <li class:active={key == sort_by}>
        <button on:click={() => (sort_by = key)}>{@html label ?? key}</button>
      </li>
    {/each}
  </ul>

  <ol>
    {#each models as [key, metadata] (key)}
      <li
        animate:flip={{ duration: 400 }}
        in:fade|local={{ delay: 100 }}
        out:fade|local={{ delay: 100 }}
      >
        <ModelCard {key} data={metadata} {stats} {sort_by} />
      </li>
    {/each}
  </ol>
</div>

<h2>Model Run Times</h2>

<p>
  Creating this benchmark (excluding debugging runs) used a total of 3137 hours of compute
  time (mix of CPU and GPU, mostly CPU). Notably, the vast majority of that was used in
  the Bayesian optimization step of the BOWSR+MEGnet model.
</p>

{#if typeof document !== `undefined`}
  <RunTimePie />
{/if}

<style>
  :is(ul, ol) {
    padding: 0;
    list-style: none;
  }
  ul {
    display: flex;
    flex-wrap: wrap;
    gap: 9pt;
    margin: 1em auto 2em;
    place-content: center;
  }
  ul > li > button {
    transition: all 0.2s;
    background-color: rgba(255, 255, 255, 0.1);
  }
  ul > li.active > button {
    background-color: darkcyan;
  }
  ol {
    display: grid;
    gap: 2em;
    grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
  }
  ol > li {
    background-color: rgba(255, 255, 255, 0.05);
    padding: 3pt 10pt 7pt;
    border-radius: 3pt;
    display: grid;
    align-content: space-between;
  }
  span {
    display: flex;
    gap: 5pt;
    place-items: center;
    place-content: center;
  }
  span :global(div.zoo-radio-btn span) {
    padding: 1pt 4pt;
  }
</style>
