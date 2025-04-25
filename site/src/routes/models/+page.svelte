<script lang="ts">
  import { ModelCard, type Metric } from '$lib'
  import { METADATA_COLS, METRICS } from '$lib/labels'
  import { get_nested_value, metric_better_as, sort_models } from '$lib/metrics'
  import { model_is_compliant, MODELS } from '$lib/models.svelte'
  import { interpolateCividis as cividis } from 'd3-scale-chromatic'
  import { ColorBar } from 'elementari'
  import { RadioButtons, Tooltip } from 'svelte-zoo'
  import { flip } from 'svelte/animate'
  import { fade } from 'svelte/transition'

  let sort_by: Metric = $state(METRICS.CPS)
  let show_non_compliant: boolean = $state(true)
  let show_details: boolean = $state(false)
  let order: `asc` | `desc` = $state(`desc`)
  let show_n_best: number = $state(MODELS.length) // show only best models
  const min_models: number = 2
  let sort_by_path: string = $derived(
    `${sort_by.path ?? ``}.${sort_by.key}`.replace(/^\./, ``),
  )

  const metric_keys = [
    `CPS`,
    `Accuracy`,
    `DAF`,
    `F1`,
    `MAE`,
    `Precision`,
    `R2`,
    `RMSE`,
    `TNR`,
    `TPR`,
    `κ_SRME`,
  ] as const
  const metrics = metric_keys.map((key) => METRICS[key])

  export const snapshot = {
    capture: () => ({ show_details, sort_by, order, show_n_best }),
    restore: (values) => ({ show_details, sort_by, order, show_n_best } = values),
  }

  function bg_color(val: number, min: number, max: number) {
    if (isNaN(val)) return `rgba(255, 255, 255, 0.1)` // Default background for NaN values
    return cividis(1 - (val - min) / (max - min)).replace(`)`, `, 0.5)`)
  }

  let lower_is_better = $derived(metric_better_as(sort_by.key) === `lower`)

  let models = $derived(
    MODELS.filter((model) => show_non_compliant || model_is_compliant(model)).sort(
      sort_models(sort_by_path, order),
    ),
  )

  let [min_val, max_val] = $derived.by(() => {
    if (sort_by.key == METADATA_COLS.model_name.key) return [] // dummy range, no color gradient needed

    // Use the helper function to get values for color scaling
    const vals = models
      .map((model) => get_nested_value(model, sort_by_path))
      .filter((val) => typeof val === `number` && !isNaN(val))
      .sort() as number[]

    return [vals.at(0) ?? 0, vals.at(-1) ?? 1]
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
    {#each [{ ...METADATA_COLS.model_name, label: `Model Name` }, ...metrics] as prop (prop.key)}
      {@const { key, label, short, description } = prop}
      <li class:active={prop.key == sort_by.key}>
        <button
          id={prop.key}
          onclick={() => {
            // Handle the case where key is 'model_name'
            sort_by = prop
            if (key === `model_name`) order = `asc`
            else order = metric_better_as(key) === `lower` ? `asc` : `desc`
          }}
          style="font-size: large;"
        >
          {@html short ?? label ?? key}
        </button>
        {#if description}
          <Tooltip text={description} max_width="20em">
            <span style="position: absolute; top: -11pt; left: -6pt; opacity: 0.6;">
              <svg><use href="#icon-info"></use></svg>
            </span>
          </Tooltip>
        {/if}
      </li>
    {/each}
  </ul>

  <legend style:opacity={sort_by.key === METADATA_COLS.model_name.key ? 0 : 1}>
    <span>
      {lower_is_better ? `best` : `worst`}
    </span>
    <ColorBar
      label="Model names colored by {sort_by.label}"
      label_style="font-size: 1.4em;"
      color_scale={cividis}
      style="min-width: min(70vw, 400px);"
      range={[min_val ?? 0, max_val ?? 1]}
    />
    <span>
      {lower_is_better ? `worst` : `best`}
    </span>
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
          {metrics}
          sort_by={sort_by_path}
          bind:show_details
          style="background-color: {bg_color(
            sort_by.key === METADATA_COLS.model_name.key
              ? 0 // No gradient for model names
              : (get_nested_value(model, sort_by_path) as number),
            lower_is_better ? max_val : min_val,
            lower_is_better ? min_val : max_val,
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
    align-items: end;
    transition: opacity 0.4s;
  }
  legend span {
    transform: translateY(7px);
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
