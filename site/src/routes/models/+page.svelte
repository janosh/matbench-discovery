<script lang="ts">
  import type { ModelStatLabel, ModelStats } from '$lib'
  import { model_is_compliant, ModelCard, MODELS } from '$lib'
  import { get_metric_value, metric_better_as } from '$lib/metrics'
  import { interpolateCividis as cividis } from 'd3-scale-chromatic'
  import { ColorBar } from 'elementari'
  import { RadioButtons, Tooltip } from 'svelte-zoo'
  import { flip } from 'svelte/animate'
  import { fade } from 'svelte/transition'
  import type { Snapshot } from './$types'

  let sort_by: keyof ModelStats | `model_name` = $state(`CPS`)
  let show_non_compliant: boolean = $state(true)
  let show_details: boolean = $state(false)
  let order: `asc` | `desc` = $state(`desc`)
  let show_n_best: number = $state(MODELS.length) // show only best models
  const min_models: number = 2

  const stats: ModelStatLabel[] = [
    { key: `CPS`, label: `CPS`, tooltip: `Combined Performance Score` },
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
    if (isNaN(val)) return `rgba(255, 255, 255, 0.1)` // Default background for NaN values
    return cividis(1 - (val - min) / (max - min)).replace(`)`, `, 0.5)`)
  }

  $effect(() => {
    // Determine default sort order based on whether lower or higher is better for this metric
    if (sort_by === `model_name`) {
      order = `asc` // Model names default to ascending alphabetical order
    } else {
      order = metric_better_as(sort_by) === `lower` ? `asc` : `desc`
    }
  })

  let sort_factor = $derived({ asc: -1, desc: 1 }[order])

  let models = $derived(
    MODELS.filter((model) => show_non_compliant || model_is_compliant(model)).sort(
      (model_1, model_2) => {
        // Special case for model_name sorting
        if (sort_by === `model_name`) {
          // For model_name, directly use localeCompare with sort_factor
          return sort_factor * model_1.model_name.localeCompare(model_2.model_name)
        }

        // Get values using the helper function for other metrics
        const val_1 = get_metric_value(model_1, sort_by)
        const val_2 = get_metric_value(model_2, sort_by)

        // Handle null, undefined, or NaN values by sorting last
        if (val_1 == null && val_2 == null) return 0
        if (val_1 == null || Number.isNaN(val_1)) return 1 // Always sort nulls/NaN to the end
        if (val_2 == null || Number.isNaN(val_2)) return -1 // Always sort nulls/NaN to the end

        if (typeof val_1 == `string` && typeof val_2 == `string`) {
          return sort_factor * (val_1 as string).localeCompare(val_2 as string)
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
      },
    ),
  )

  let [min_val, max_val] = $derived.by(() => {
    if (sort_by === `model_name`) return [0, 1] // dummy range, no color gradient needed

    // Use the helper function to get values for color scaling
    const vals = models
      .map((model) => get_metric_value(model, sort_by, `unique_prototypes`))
      .filter((val) => typeof val === `number` && !isNaN(val))
      .sort() as number[]

    // Determine color range based on whether lower or higher is better
    const lower_better = metric_better_as(sort_by) === `lower`
    const [min, max] = [vals.at(0), vals.at(-1)]
    return lower_better ? [max, min] : [min, max]
  })
</script>

<div style="display: grid;">
  <span>
    <input type="checkbox" bind:checked={show_non_compliant} />Show non-compliant models
    &ensp; &emsp;&emsp; Sort
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
            // Handle the case where key is 'model_name'
            if (key === `model_name`) {
              sort_by = `model_name`
              order = `asc` // Default for model names
            } else {
              sort_by = key as keyof ModelStats
              order = metric_better_as(key) === `lower` ? `asc` : `desc`
            }
          }}
          style="font-size: large;"
        >
          {@html label ?? key}
        </button>
        {#if tooltip}
          <Tooltip text={tooltip} max_width="20em">
            <span style="position: absolute; top: -11pt; left: -6pt; opacity: 0.6;">
              <svg><use href="#icon-info"></use></svg>
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
        style="grid-row: span {show_details ? 5 : 4};"
      >
        <ModelCard
          {model}
          {stats}
          {sort_by}
          bind:show_details
          style="background-color: {bg_color(
            sort_by === `model_name`
              ? 0 // No gradient for model names
              : (get_metric_value(model, sort_by, `unique_prototypes`) ?? 0),
            min_val ?? 0,
            max_val ?? 1,
          )};"
        />
      </li>
    {/each}
  </ol>
</div>

<!-- link to ALL model pages with hidden links for the SvelteKit crawler -->
{#each MODELS as model (model.model_name)}
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
    gap: 1em;
    grid-template-columns: repeat(auto-fit, minmax(420px, 1fr));
  }
  ol > li {
    background-color: rgba(255, 255, 255, 0.05);
    padding: 6pt 10pt 14pt;
    border-radius: 3pt;
    display: grid;
    grid-template-rows: subgrid;
    position: relative;
    gap: 1em;
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
    transform: scale(1.2);
    padding: 2pt;
  }
  input[type='number']::-webkit-inner-spin-button {
    display: none;
  }
  input[type='checkbox'] {
    transform: scale(1.3);
  }
</style>
