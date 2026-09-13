<script lang="ts">
  import ModelSelect from '$lib/ModelSelect.svelte'
  import { ACTIVE_MODELS } from '$lib/models.svelte'
  import { max } from 'd3-array'
  import { ColorBar } from 'matterviz/plot'
  import { format_num } from 'matterviz/labels'
  import { PeriodicTable, TableInset } from 'matterviz/periodic-table'
  import { per_element_each_errors as each_errors } from '$lib/per-element-errors'
  import { bind_url_params } from '$lib/url-state.svelte'
  import { bool_from_param, bool_url_entry } from 'svelte-widgets/url-params'

  const model_options = ACTIVE_MODELS.filter(
    ({ model_key }) => model_key in each_errors,
  ).map(({ model_key, model_name }) => ({ label: model_name, value: model_key }))
  let selected_models = $state(model_options.slice(0, 1))
  let manual_cbar_max = $state(false)
  let normalized = $state(true)
  let cbar_max = $state(0.3)

  // where each model's value lands in a split tile, keyed by selection count
  // (matches matterviz ElementTile's auto layouts: 2=diagonal, 3=horizontal, 4=quadrant)
  const split_positions: Record<number, string[]> = {
    2: [`top left`, `bottom right`],
    3: [`top`, `middle`, `bottom`],
    4: [`top left`, `top right`, `bottom left`, `bottom right`],
  }

  bind_url_params(
    (params) => {
      const keys = new Set(params.get(`element_models`)?.split(`,`))
      selected_models = [...keys]
        .flatMap((key) => model_options.find(({ value }) => value === key) ?? [])
        .slice(0, 4)
      if (!selected_models.length) selected_models = model_options.slice(0, 1)
      normalized = bool_from_param(params, `element_normalized`, true)
      manual_cbar_max = bool_from_param(params, `element_manual_max`)
      const maximum = Number(params.get(`element_max`))
      cbar_max = maximum >= 0.01 && maximum <= 0.7 ? maximum : 0.3
    },
    () => [
      [
        `element_models`,
        selected_models.map(({ value }) => value).join(`,`),
        model_options[0]?.value ?? ``,
      ],
      bool_url_entry(`element_normalized`, normalized, true),
      bool_url_entry(`element_manual_max`, manual_cbar_max),
      [`element_max`, String(cbar_max), `0.3`],
    ],
  )

  const test_set_std = each_errors[`Test set standard deviation`]

  // Selection order determines segment order; one value paints a solid tile.
  let heatmap_values = $derived(
    Object.fromEntries(
      Object.entries(test_set_std).map(([element, std]) => [
        element,
        selected_models.map(({ value }) => {
          const error = each_errors[value][element]
          const denom = normalized ? std : 1
          return error == null || !denom ? `n/a` : error / denom
        }),
      ]),
    ),
  )
  // Non-numeric labels use the table's missing style and never enter the color scale.
  let current_data_max = $derived(
    max(Object.values(heatmap_values).flat(), (value) =>
      typeof value === `number` ? value : undefined,
    ) ?? 0,
  )
  let cs_range = $derived<[number, number]>([
    0,
    manual_cbar_max ? cbar_max : current_data_max,
  ])
  let error_unit = $derived(normalized ? `normalized` : `eV/atom`)
  let cbar_title = $derived(
    `${selected_models.length === 1 ? selected_models[0].label : `Element-projected error`} (${error_unit})`,
  )
</script>

<p>
  This periodic table heatmap shows the MAE of model-predicted convex hull distance
  projected onto each element. The errors for every structure in the test set are
  projected onto the fraction of each element in the composition and averaged over all
  structures. The error is the absolute difference per atom between predicted and actual
  energy distance to the convex hull. Select up to 4 models to compare them side by side:
  each element tile splits into one segment per model.
</p>

<ModelSelect
  bind:value={selected_models}
  options={model_options}
  max_select={4}
  min_select={1}
/>

{#if selected_models.length > 1}
  <div class="split-legend">
    {#each selected_models as model, idx (model.value)}
      <span>
        <strong>{model.label}</strong>
        <small>({split_positions[selected_models.length]?.[idx]})</small>
      </span>
    {/each}
  </div>
{/if}

<form>
  <label>
    Manual color bar max
    <input type="checkbox" bind:checked={manual_cbar_max} />
    <input
      type="range"
      disabled={!manual_cbar_max}
      bind:value={cbar_max}
      min={0.01}
      max={0.7}
      step={0.001}
    />
    {cbar_max}
  </label>
  <label>
    <input type="checkbox" bind:checked={normalized} />
    Divide each element value by its std. dev. of target energies over all test structures containing
    a given element
  </label>
  <small>
    This is meant to correct for the fact that some elements are inherently more difficult
    to predict since some have a wider distribution of convex hull distances.
  </small>
</form>

<PeriodicTable
  {heatmap_values}
  color_scale="interpolateViridis"
  color_scale_range={cs_range}
  tile_props={{ float_fmt: `.2` }}
  show_photo={false}
  missing={{ color: `rgba(255,255,255,0.3)` }}
>
  {#snippet inset({ active_element })}
    {#if active_element}
      <TableInset style="align-content: center">
        <strong class="model-errors">
          {active_element.name}: {#each selected_models as model, idx (model.value)}
            {#if idx > 0}&ensp;{/if}
            {@const elem_error = heatmap_values[active_element.symbol]?.[idx] ?? `n/a`}
            <span>
              {#if selected_models.length > 1}{model.label}{/if}
              <b>{typeof elem_error === `number` ? format_num(elem_error) : elem_error}</b
              >
            </span>
          {/each}
          <small>{error_unit}</small>
        </strong>
        <ColorBar
          title={cbar_title}
          title_side="top"
          scale="interpolateViridis"
          tick_labels={5}
          range={cs_range}
          style="width: 85%; margin: 0 2em"
        />
      </TableInset>
    {/if}
  {/snippet}
</PeriodicTable>

<style>
  form {
    display: flex;
    flex-direction: column;
    margin: 1em;
    gap: 10pt;
  }
  form label {
    display: flex;
    place-content: center;
    gap: 1ex;
  }
  form label + small {
    max-width: 60em;
    margin: 0 auto;
    text-align: center;
  }
  .split-legend {
    display: flex;
    flex-wrap: wrap;
    place-content: center;
    gap: 1em;
    margin: 1ex auto;
  }
  .split-legend small {
    color: var(--text-secondary);
  }
  .model-errors {
    display: block;
    text-align: center;
    min-height: 18pt;
  }
  .model-errors span,
  .model-errors small {
    font-weight: lighter;
  }
</style>
