<script lang="ts">
  import type { ModelStatLabel, ModelStats } from '$lib'
  import { model_is_compliant, MODEL_METADATA, ModelCard } from '$lib'
  import { discovery } from '$pkg/modeling-tasks.yml'
  import { interpolateCividis as cividis } from 'd3-scale-chromatic'
  import { ColorBar } from 'elementari'
  import 'iconify-icon'
  import { RadioButtons, Toggle, Tooltip } from 'svelte-zoo'
  import { flip } from 'svelte/animate'
  import { fade } from 'svelte/transition'
  import type { Snapshot } from './$types'

  let sort_by: keyof ModelStats | `model_name` = $state(`F1`)
  let show_non_compliant: boolean = $state(true)
  let show_details: boolean = $state(false)
  let order: `asc` | `desc` = $state(`desc`)
  let show_n_best: number = $state(MODEL_METADATA.length) // show only best models
  const min_models: number = 2

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
    {
      key: `κ_SRME`,
      label: `κ<sub>SRME</sub>`,
      tooltip: `symmetric relative mean error in predicted phonon mode contributions to thermal conductivity κ`,
    },
  ]

  export const snapshot: Snapshot = {
    capture: () => ({ show_details, sort_by, order, show_n_best }),
    restore: (values) => ({ show_details, sort_by, order, show_n_best } = values),
  }

  function bg_color(val: number, min: number, max: number) {
    return cividis(1 - (val - min) / (max - min)).replace(`)`, `, 0.5)`)
  }
  $effect(() => {
    order = discovery.metrics.lower_is_better.includes(sort_by) ? `asc` : `desc`
  })
  let sort_factor = $derived({ asc: -1, desc: 1 }[order])
  let models = $derived(
    MODEL_METADATA.filter(
      (model) => show_non_compliant || model_is_compliant(model),
    ).sort((model_1, model_2) => {
      const metrics_1 = model_1.metrics?.discovery?.full_test_set ?? {}
      const metrics_2 = model_2.metrics?.discovery?.full_test_set ?? {}
      const [val_1, val_2] = [metrics_1[sort_by], metrics_2[sort_by]]

      // Handle null values by sorting last
      if (val_1 == null && val_2 == null) return 0
      if (val_1 == null) return 1
      if (val_2 == null) return -1

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
    }),
  )
  let [min_val, max_val] = $derived.by(() => {
    const vals = models
      .map((model) => model.metrics?.discovery?.full_test_set?.[sort_by])
      .filter((val) => typeof val === `number`)

    const lower_better = discovery.metrics.lower_is_better.includes(sort_by)
    return lower_better
      ? [Math.max(...vals), Math.min(...vals)]
      : [Math.min(...vals), Math.max(...vals)]
  })
</script>

<div style="display: grid;">
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
    {#each [{ key: `model_name`, label: `Model Name`, tooltip: undefined }, ...stats] as { key, label, tooltip } (key)}
      <li class:active={key == sort_by}>
        <button
          id={key}
          onclick={() => {
            sort_by = key
            order = discovery.metrics.lower_is_better.includes(key) ? `asc` : `desc`
          }}
          style="font-size: large;"
        >
          {@html label ?? key}
        </button>
        {#if tooltip}
          <Tooltip text={tooltip} max_width="20em">
            <span style="position: absolute; top: -1ex; left: -4pt; color: gray;">
              <iconify-icon icon="octicon:info-16" aria-label="Info" height="9.5pt"
              ></iconify-icon>
            </span>
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

  <ol class="models">
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
            model.metrics?.discovery?.unique_prototypes?.[sort_by],
            min_val,
            max_val,
          )};"
        />
      </li>
    {/each}
  </ol>
</div>

<!-- link to ALL model pages with hidden links for the SvelteKit crawler -->
{#each MODEL_METADATA as model (model.model_name)}
  <a href="/models/{model.model_key}" hidden>
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
