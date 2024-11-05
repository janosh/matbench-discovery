<script lang="ts">
  import type { ModelStatLabel, ModelStats } from '$lib'
  import { model_is_compliant, MODEL_METADATA, ModelCard } from '$lib'
  import { lower_is_better } from '$root/scripts/metrics-which-is-better.yml'
  import Icon from '@iconify/svelte'
  import { interpolateCividis as cividis } from 'd3-scale-chromatic'
  import { ColorBar } from 'elementari'
  import { RadioButtons, Toggle, Tooltip } from 'svelte-zoo'
  import { flip } from 'svelte/animate'
  import { fade } from 'svelte/transition'
  import type { Snapshot } from './$types'

  let sort_by: keyof ModelStats | `model_name` = `F1`
  let show_non_compliant: boolean = false
  let show_details: boolean = false
  let order: `asc` | `desc` = `desc`
  let show_n_best: number = MODEL_METADATA.length // show only best models
  const min_models: number = 2

  $: models = MODEL_METADATA.filter(
    (model) => show_non_compliant || model_is_compliant(model),
  ).sort((model_1, model_2) => {
    const metrics_1 = model_1.metrics?.discovery?.full_test_set ?? {}
    const metrics_2 = model_2.metrics?.discovery?.full_test_set ?? {}
    const [val_1, val_2] = [metrics_1[sort_by], metrics_2[sort_by]]

    // Handle null values by sorting last
    if (val_1 === null && val_2 === null) return 0
    if (val_1 === null) return 1
    if (val_2 === null) return -1

    if (typeof val_1 == `string`) {
      return sort_factor * val_1.localeCompare(val_2)
    } else if (typeof val_1 == `number` && typeof val_2 == `number`) {
      // interpret runt_time==0 as infinity
      if (sort_by == `Run Time (h)`) {
        if (val_1 == 0) return -sort_factor
        if (val_2 == 0) return sort_factor
      }
      return sort_factor * (val_2 - val_1)
    } else {
      throw `Unexpected type '${val_1}' encountered sorting by key '${sort_by}'`
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
    {
      key: `κ<sub>SRME</sub>`,
      tooltip: `symmetric relative mean error in predicted phonon mode contributions to thermal conductivity κ`,
    },
  ]

  export const snapshot: Snapshot = {
    capture: () => ({ show_details, sort_by, order, show_n_best }),
    restore: (values) => ({ show_details, sort_by, order, show_n_best } = values),
  }

  $: sort_factor = { asc: -1, desc: 1 }[order]
  $: min_val = Math.min(
    ...models.map((model) => model.metrics?.discovery?.full_test_set[sort_by] as number),
  )
  $: max_val = Math.max(
    ...models.map((model) => model.metrics?.discovery?.full_test_set[sort_by] as number),
  )
  $: if (lower_is_better.includes(sort_by)) [min_val, max_val] = [max_val, min_val]
  $: order = lower_is_better.includes(sort_by) ? `asc` : `desc`

  function bg_color(val: number, min: number, max: number) {
    return cividis(1 - (val - min) / (max - min)).replace(`)`, `, 0.5)`)
  }
</script>

<div style="display: grid; margin: 1vw 3vw;">
  <h1>Leaderboard</h1>

  <p style="text-align: center;">
    Sort models by different metrics (thermodynamic stability classification, convex hull
    distance regressions or tun time).
  </p>

  <span>
    <Toggle bind:checked={show_non_compliant}>Show non-compliant models&ensp;</Toggle>
    &emsp;&emsp; Sort
    <input type="number" min={min_models} max={models.length} bind:value={show_n_best} />
    best models
    <RadioButtons bind:selected={order} options={[`asc`, `desc`]} /> by:
  </span>
  <ul>
    {#each [{ key: `model_name`, label: `Model Name` }, ...stats] as { key, label, tooltip }}
      <li class:active={key == sort_by}>
        <button
          id={key}
          on:click={() => {
            sort_by = key
            order = lower_is_better.includes(key) ? `asc` : `desc`
          }}
          style="font-size: large;"
        >
          {@html label ?? key}
        </button>
        {#if tooltip}
          <Tooltip
            text={tooltip}
            tip_style="white-space: nowrap; font-size: 9pt;"
            max_width="20em"
            style="position: absolute; transform: translate(-45%, -45%); color: gray;"
          >
            <Icon icon="octicon:info-16" title="Info" height="9pt" />
          </Tooltip>
        {/if}
      </li>
    {/each}
  </ul>

  <legend>
    heading color best
    <ColorBar color_scale={cividis} style="min-width: min(70vw, 400px);" />
    worst
  </legend>

  <ol>
    {#each models.slice(0, Math.max(min_models, show_n_best)) as model (model.model_name)}
      <li
        animate:flip={{ duration: 400 }}
        in:fade={{ delay: 100 }}
        out:fade={{ delay: 100 }}
      >
        <ModelCard
          {model}
          {stats}
          {sort_by}
          bind:show_details
          style="background-color: {bg_color(
            model.metrics?.discovery?.full_test_set[sort_by],
            min_val,
            max_val,
          )};"
        />
        <!-- maybe show this text in a tooltip: This model was not trained on the
        canonical training set. It's results should not be seen as a one-to-one
        comparison to the other models but rather proof of concept of what is possible. -->
        <!-- {#if model.training_set}
          <strong class="train-set">
            <Icon icon="ion:ios-warning" inline />
            Custom training set: {model.training_set.title}
          </strong>
        {/if} -->
      </li>
    {/each}
  </ol>
</div>

<!-- link to ALL model pages with hidden links for the crawler -->
{#each MODEL_METADATA as model}
  <a href="/models/{model.model_name.toLowerCase().replaceAll(` `, `-`)}" hidden>
    {model.model_name}
  </a>
{/each}

<style>
  legend {
    display: flex;
    gap: 8pt;
    place-content: center;
    opacity: 0.8;
    margin: 2em auto;
    font-weight: lighter;
  }
  :is(ul, ol) {
    padding: 0;
    list-style: none;
  }
  ul {
    display: flex;
    flex-wrap: wrap;
    gap: 9pt;
    margin: 2.5ex auto 3ex;
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
  /* strong.train-set {
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
  } */
</style>
