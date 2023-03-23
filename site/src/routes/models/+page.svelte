<script lang="ts">
  import type { ModelStatLabel, ModelStats } from '$lib'
  import { ModelCard } from '$lib'
  import Icon from '@iconify/svelte'
  import { RadioButtons, Tooltip } from 'svelte-zoo'
  import { flip } from 'svelte/animate'
  import { fade } from 'svelte/transition'
  import type { Snapshot } from './$types'

  export let data

  let sort_by: keyof ModelStats | 'model_name' = `F1`
  let show_details: boolean = false
  let order: 'asc' | 'desc' = `desc`
  let n_best: number = 8 // show only best models
  const min_models: number = 2
  $: sort_factor = { asc: -1, desc: 1 }[order]

  $: models = data.models.sort((md1, md2) => {
    if (typeof md1[sort_by] == `string`) {
      return sort_factor * -md1[sort_by].localeCompare(md2[sort_by])
    } else if (typeof md1[sort_by] == `number`) {
      return sort_factor * (md2[sort_by] - md1[sort_by])
    } else {
      console.error(
        `Sorting by key ${sort_by} gives unknown type: ${typeof md1[sort_by]}`
      )
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
    capture: () => ({ show_details, sort_by, order }),
    restore: (values) => ({ show_details, sort_by, order } = values),
  }
</script>

<div class="pull-left">
  <h1>Models</h1>

  <span>
    Sort <input type="number" min={min_models} max={models.length} bind:value={n_best} />
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
    {#each models.slice(0, Math.max(min_models, n_best)) as data, idx (data.model_name)}
      <li
        animate:flip={{ duration: 400 }}
        in:fade|local={{ delay: 100 }}
        out:fade|local={{ delay: 100 }}
      >
        <ModelCard {idx} {data} {stats} {sort_by} bind:show_details />
      </li>
    {/each}
  </ol>
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
    grid-template-columns: repeat(auto-fit, minmax(340px, 1fr));
  }
  ol > li {
    background-color: rgba(255, 255, 255, 0.05);
    padding: 3pt 10pt 7pt;
    border-radius: 3pt;
    display: grid;
    align-content: start;
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
</style>
