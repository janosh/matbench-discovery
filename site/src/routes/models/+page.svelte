<script lang="ts">
  import { Icon, type Label, ModelCard } from '$lib'
  import { ALL_METRICS, METADATA_COLS } from '$lib/labels'
  import { get_nested_value, metric_better_as, sort_models } from '$lib/metrics'
  import { model_is_compliant, MODELS } from '$lib/models.svelte'
  import { interpolateRdBu } from 'd3-scale-chromatic'
  import { ColorBar, luminance } from 'matterviz'
  import { RadioButtons } from 'svelte-multiselect'
  import { tooltip } from 'svelte-multiselect/attachments'
  import { flip } from 'svelte/animate'
  import { fade } from 'svelte/transition'

  let sort_by: Label = $state(ALL_METRICS.CPS)
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
  const metrics = metric_keys.map((key) => ALL_METRICS[key])

  export const snapshot = {
    capture: () => ({ show_details, sort_by, order, show_n_best }),
    restore: (values) => ({ show_details, sort_by, order, show_n_best } = values),
  }

  function bg_color(val: number, min: number, max: number) {
    if (isNaN(val)) return `rgba(255, 255, 255, 0.6)` // Default background for NaN values
    return interpolateRdBu((val - min) / (max - min))
  }

  let lower_is_better = $derived(metric_better_as(sort_by.key) === `lower`)

  let models = $derived(
    MODELS.filter((model) => show_non_compliant || model_is_compliant(model)).sort(
      sort_models(sort_by_path, order),
    ),
  )

  let [min_val, max_val] = $derived.by(() => {
    if (!sort_by.better) return [NaN, NaN]

    const vals = models
      .map((model) => get_nested_value(model, sort_by_path))
      .filter((val) => typeof val === `number` && !isNaN(val))
      .sort() as number[]

    return [vals.at(0) ?? 0, vals.at(-1) ?? 1]
  })
  let [best_val, worst_val] = $derived(
    lower_is_better ? [max_val, min_val] : [min_val, max_val],
  )
</script>

<div style="display: grid">
  <span>
    <input type="checkbox" bind:checked={show_non_compliant} />Show non-compliant models
    &ensp; &emsp;&emsp; Sort
    <input type="number" min={min_models} max={models.length} bind:value={show_n_best} />
    best models
    <RadioButtons bind:selected={order} options={[`asc`, `desc`]} /> by:
  </span>

  <ul>
    {#each [{ ...METADATA_COLS.model_name, label: `Model Name` }, ...metrics] as
      prop
      (prop.key)
    }
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
          style="font-size: large; height: 26pt; position: relative"
        >
          {@html short ?? label ?? key}
          {#if description}
            <span
              {@attach tooltip({ content: description })}
              style="width: 14pt; height: 14pt; position: absolute; top: -7pt; right: -7pt; opacity: 0.6"
            >
              <Icon icon="Info" />
            </span>
          {/if}
        </button>
      </li>
    {/each}
  </ul>

  <legend style:opacity={sort_by.key === METADATA_COLS.model_name.key ? 0 : 1}>
    <span>
      {lower_is_better ? `best` : `worst`}
    </span>
    <ColorBar
      title="Card titles colored by {sort_by.label}"
      title_style="font-size: 1.5em;"
      color_scale={lower_is_better ? (t) => interpolateRdBu(1 - t) : interpolateRdBu}
      style="min-width: min(70vw, 400px); height: 14pt"
      range={lower_is_better ? [worst_val, best_val] : [best_val, worst_val]}
    />
    <span>
      {lower_is_better ? `worst` : `best`}
    </span>
  </legend>

  <ol class="models full-bleed">
    {#each models.slice(0, Math.max(min_models, show_n_best)) as model (model.model_name)}
      {@const metric_val = sort_by.better ? get_nested_value(model, sort_by_path) : 0}
      {@const bg_clr = bg_color(metric_val as number, best_val, worst_val)}
      {@const text_color = luminance(bg_clr) > 0.7 ? `black` : `white`}
      <li
        animate:flip={{ duration: 400 }}
        in:fade={{ delay: 100 }}
        out:fade={{ delay: 100 }}
        style:grid-row="span {show_details ? 5 : 4}"
      >
        <ModelCard
          {model}
          {metrics}
          sort_by={sort_by_path}
          bind:show_details
          title_style="background-color: {bg_clr}; color: {text_color};"
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
    transform: translateY(4px);
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
  }
  ul > li.active button {
    background-color: var(--btn-bg);
  }
  ol {
    display: grid;
    gap: 1em;
    grid-template-columns: repeat(auto-fit, minmax(420px, 1fr));
  }
  ol > li {
    background-color: var(--blockquote-bg);
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
