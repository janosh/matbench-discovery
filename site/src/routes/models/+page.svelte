<script lang="ts">
  import type { ModelStatLabel, ModelStats } from '$lib'
  import { ModelCard } from '$lib'
  import Icon from '@iconify/svelte'
  import { RadioButtons, Tooltip } from 'svelte-zoo'
  import { flip } from 'svelte/animate'
  import { fade } from 'svelte/transition'
  import type { Snapshot } from './$types'
  import ElementErrorsPtableHeatmap from './element-errors-ptable-heatmap.svelte'

  export let data

  let sort_by: keyof ModelStats | 'model_name' = `F1`
  let show_details: boolean = false
  let order: 'asc' | 'desc' = `desc`
  let show_n_best: number = 8 // show only best models
  const min_models: number = 2
  $: sort_factor = { asc: -1, desc: 1 }[order]

  $: models = data.models.sort((model_1, model_2) => {
    const [val_1, val_2] = [model_1[sort_by], model_2[sort_by]]
    if (typeof val_1 == `string`) {
      return sort_factor * val_1.localeCompare(val_2)
    } else if (typeof val_1 == `number`) {
      return sort_factor * (val_2 - val_1)
    } else {
      console.error(`Sorting by key ${sort_by} gives unknown type: ${typeof val_1}`)
    }
  })
  const stats: ModelStatLabel[] = [
    { key: `Accuracy` },
    { key: `DAF`, tooltip: `Discovery Acceleration Factor` },
    { key: `F1` },
    { key: `MAE`, unit: `eV / atom`, tooltip: `Mean Absolute Error` },
    { key: `Precision` },
    { key: `R2`, label: `R<sup>2</sup>` },
    { key: `RMSE`, unit: `eV / atom`, tooltip: `Root Mean Squared Error` },
    { key: `TNR`, tooltip: `True Negative Rate` },
    { key: `TPR`, tooltip: `True Positive Rate` },
    { key: `Run Time (h)`, label: `Run time`, unit: `h` },
  ]

  export const snapshot: Snapshot = {
    capture: () => ({ show_details, sort_by, order, show_n_best }),
    restore: (values) => ({ show_details, sort_by, order, show_n_best } = values),
  }
</script>

<div>
  <h1>Models</h1>

  <span>
    Sort
    <input type="number" min={min_models} max={models.length} bind:value={show_n_best} />
    best models
    <RadioButtons bind:selected={order} options={[`asc`, `desc`]} /> by:
  </span>
  <ul>
    {#each [{ key: `model_name`, label: `Model Name` }, ...stats] as { key, label, tooltip }}
      <li class:active={key == sort_by}>
        <button id={key} on:click={() => (sort_by = key)}>{@html label ?? key}</button>
        {#if tooltip}
          <Tooltip
            text={tooltip}
            tip_style="white-space: nowrap; font-size: 9pt;"
            max_width="20em"
            style="position: absolute; transform: translate(-45%, -45%); color: gray;"
          >
            <Icon icon="material-symbols:info-outline" title="Info" height="9pt" />
          </Tooltip>
        {/if}
      </li>
    {/each}
  </ul>

  <ol>
    {#each models.slice(0, Math.max(min_models, show_n_best)) as data (data.model_name)}
      <li
        animate:flip={{ duration: 400 }}
        in:fade|local={{ delay: 100 }}
        out:fade|local={{ delay: 100 }}
      >
        <ModelCard {data} {stats} {sort_by} bind:show_details />
        {#if data.training_set}
          <!-- maybe show this text in a tooltip: This model was not trained on the
            canonical training set. It's results should not be seen as a one-to-one
            comparison to the other models but rather proof of concept of what is possible. -->
          <strong class="train-set">
            <Icon icon="ion:ios-warning" inline />
            Custom training set: {data.training_set}
          </strong>
        {/if}
      </li>
    {/each}
  </ol>
</div>

<div style="margin-top: 6em;">
  <h2>Per-Element Model Error Heatmaps</h2>

  <ElementErrorsPtableHeatmap />
</div>

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
  ul > li button {
    transition: all 0.2s;
    background-color: rgba(255, 255, 255, 0.1);
  }
  ul > li.active button {
    background-color: darkcyan;
  }
  ol {
    display: grid;
    gap: 2em;
    grid-template-columns: repeat(auto-fit, minmax(360px, 1fr));
  }
  ol > li {
    background-color: rgba(255, 255, 255, 0.05);
    padding: 6pt 10pt 14pt;
    border-radius: 3pt;
    display: grid;
    align-content: start;
    position: relative;
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
  input[type='number'] {
    text-align: center;
  }
  input[type='number']::-webkit-inner-spin-button {
    display: none;
  }
  strong.train-set {
    display: block;
    background-color: rgb(174, 79, 28);
    color: white;
    position: absolute;
    left: 50%;
    transform: translate(-50%, 50%);
    bottom: 0;
    padding: 1pt 3pt;
    border-radius: 1ex;
    font-size: smaller;
  }
</style>
