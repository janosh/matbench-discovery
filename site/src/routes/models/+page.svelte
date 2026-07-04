<script lang="ts">
  import { type Label, ModelCard } from '$lib'
  import { ALL_METRICS, MD_METRICS, METADATA_COLS } from '$lib/labels'
  import { get_nested_value, metric_better_as, sort_models } from '$lib/metrics'
  import { model_is_compliant, MODELS } from '$lib/models.svelte'
  import { bind_url_params, valid_query_param } from '$lib/url-state.svelte'
  import { interpolateRdBu } from 'd3-scale-chromatic'
  import { ColorBar, Icon, pick_contrast_color } from 'matterviz'
  import { untrack } from 'svelte'
  import { tooltip } from 'svelte-multiselect/attachments'
  import { flip } from 'svelte/animate'
  import { fade } from 'svelte/transition'

  // Accept data prop for SvelteKit compliance (used for testing initial_show_n_best)
  let { data }: { data?: { initial_show_n_best?: number } } = $props()

  let sort_by: Label = $state(ALL_METRICS.CPS)
  let show_non_compliant: boolean = $state(true)
  let show_details: boolean = $state(false)
  let order: `asc` | `desc` = $state(`desc`)
  const min_models: number = 2
  // Enforce minimum and maximum when initializing from prop (intentionally captures initial value only)
  const initial_n_best = untrack(() =>
    Math.min(
      MODELS.length,
      Math.max(min_models, data?.initial_show_n_best ?? MODELS.length),
    ),
  )
  let show_n_best: number = $state(initial_n_best)
  let sort_by_path: string = $derived(
    `${sort_by.path ?? ``}.${sort_by.key}`.replace(/^\./, ``),
  )

  // only the 1-2 headline metrics per task (CPS overall; F1+MAE discovery; RMSD
  // geo-opt; κ_SRME phonons; CMDS+ΔvDOS molecular dynamics) -- the full metric set
  // lives in the landing-page table, this sort list stays skimmable
  const metric_keys = [`CPS`, `F1`, `MAE`, `RMSD`, `κ_SRME`] as const
  const metrics = [
    ...metric_keys.map((key) => ALL_METRICS[key]),
    MD_METRICS.md_combined_score,
    MD_METRICS.md_vdos_error,
  ]
  const sort_options = [{ ...METADATA_COLS.model_name, label: `Model Name` }, ...metrics]

  const sort_keys = new Set(sort_options.map((opt) => opt.key))
  const order_dirs = new Set([`asc`, `desc`] as const)
  const default_sort_key = ALL_METRICS.CPS.key
  const read_url_params = (params: URLSearchParams) => {
    const sort_key = valid_query_param(params, `sort`, default_sort_key, sort_keys)
    sort_by = sort_options.find((opt) => opt.key === sort_key) ?? ALL_METRICS.CPS
    order = valid_query_param(params, `dir`, `desc`, order_dirs)
    // any non-numeric or out-of-range n_best (incl. missing param -> 0) resets
    const n_best = Math.trunc(Number(params.get(`n_best`) ?? ``))
    show_n_best = n_best >= min_models ? Math.min(n_best, MODELS.length) : initial_n_best
    show_non_compliant = params.get(`non_compliant`) !== `0`
  }
  bind_url_params(read_url_params, () => [
    [`sort`, sort_by.key, default_sort_key],
    [`dir`, order, `desc`],
    [`n_best`, `${show_n_best}`, `${initial_n_best}`],
    [`non_compliant`, show_non_compliant ? `` : `0`],
  ])

  const capture_state = () => ({ show_details, sort_by, order, show_n_best })
  export const snapshot = {
    capture: capture_state,
    restore: (values: ReturnType<typeof capture_state>) =>
      ({ show_details, sort_by, order, show_n_best } = values),
  }

  function bg_color(val: number, min: number, max: number) {
    if (isNaN(val)) return `rgba(255, 255, 255, 0.6)` // Default background for NaN values
    return interpolateRdBu((val - min) / (max - min))
  }

  let lower_is_better = $derived(metric_better_as(sort_by.key) === `lower`)

  let models = $derived(
    MODELS.filter((model) => show_non_compliant || model_is_compliant(model)).toSorted(
      sort_models(sort_by_path, order),
    ),
  )

  let [min_val, max_val] = $derived.by(() => {
    if (!sort_by.better) return [NaN, NaN]

    const vals = models
      .map((model) => get_nested_value(model, sort_by_path))
      .filter((val): val is number => typeof val === `number` && !isNaN(val))

    return vals.length ? [Math.min(...vals), Math.max(...vals)] : [0, 1]
  })
  let [best_val, worst_val] = $derived(
    lower_is_better ? [max_val, min_val] : [min_val, max_val],
  )
</script>

<!-- explicit minmax(0, 1fr) column: with the default auto track, the full-bleed
ol's own 100vw-based width inflates the track beyond this wrapper's content width
(% margins count as 0 during track sizing), so the ol's -50vw + 50% centering
resolved against the wrong containing block -> off-center at 100% zoom and
horizontal overflow at wider viewports (e.g. browser zoom < 100%) -->
<div style="display: grid; grid-template-columns: minmax(0, 1fr)">
  <span>
    <input type="checkbox" bind:checked={show_non_compliant} />Show non-compliant models
    &ensp; &emsp;&emsp; Sort
    <input type="number" min={min_models} max={models.length} bind:value={show_n_best} />
    best models
    <span class="radio-group">
      {#each [`asc`, `desc`] as value (value)}
        <label>
          <input type="radio" name="order" {value} bind:group={order} />
          {value}
        </label>
      {/each}
    </span> by:
  </span>

  <ul>
    {#each sort_options as prop (prop.key)}
      {@const { key, label, description } = prop}
      <li class:active={prop.key == sort_by.key}>
        <button
          id={prop.key}
          onclick={() => {
            sort_by = prop
            // ascending for model name (alphabetical) and lower=better metrics
            order = key === `Model` || metric_better_as(key) === `lower` ? `asc` : `desc`
          }}
          style="position: relative"
        >
          {@html label ?? key}
          {#if description}
            <!-- round page-bg backing so the glyph's transparent cutouts don't show
            the button's corner behind the badge (no opacity for the same reason) -->
            <span
              {@attach tooltip({ content: description, allow_html: true })}
              style="width: 8pt; height: 8pt; position: absolute; top: -4pt; right: -4pt; background: var(--page-bg); border-radius: 50%"
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
      style="min-width: min(70vw, 400px)"
      bar_style="height: 14pt;"
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
      {@const text_color = pick_contrast_color({ bg_color: bg_clr })}
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
  /* plain --btn-bg was invisible here since it's already every button's default
  background; blue tint + inset ring (theme link color) marks the active sort */
  ul > li.active button {
    background-color: color-mix(in srgb, var(--link-color) 25%, transparent);
    box-shadow: inset 0 0 0 1px var(--link-color);
  }
  ol {
    display: grid;
    gap: 1em;
    grid-template-columns: repeat(auto-fit, minmax(420px, 1fr));
  }
  ol > li {
    background-color: light-dark(var(--light-surface-bg), var(--dark-blockquote-bg));
    border: 1px solid var(--card-border);
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
  .radio-group {
    gap: 5pt;
  }
  input[type='number'] {
    text-align: center;
    padding: 2pt;
  }
  input[type='number']::-webkit-inner-spin-button {
    display: none;
  }
  input[type='checkbox'] {
    transform: scale(1.3);
  }
</style>
