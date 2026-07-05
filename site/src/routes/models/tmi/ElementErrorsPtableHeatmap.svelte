<script lang="ts">
  import { MODELS, PtableInset } from '$lib'
  import type { ModelData } from '$lib/types'
  import type { ChemicalElement, ElementSymbol } from 'matterviz'
  import { ColorBar, format_num, PeriodicTable, TableInset } from 'matterviz'
  import type { D3InterpolateName } from 'matterviz/colors'
  import type { ComponentProps } from 'svelte'
  import Select from 'svelte-multiselect'
  import { per_element_each_errors as each_errors } from '$lib/per-element-errors'

  const models_with_errors = MODELS.filter(
    (model): model is ModelData & { model_key: string } =>
      typeof model.model_key === `string` && model.model_key in each_errors,
  )
  const model_names_with_errors = models_with_errors.map((model) => model.model_name)

  // where each model's value lands in a split tile, keyed by selection count
  // (matches matterviz ElementTile's auto layouts: 2=diagonal, 3=horizontal, 4=quadrant)
  const split_positions: Record<number, string[]> = {
    2: [`top left`, `bottom right`],
    3: [`top`, `middle`, `bottom`],
    4: [`top left`, `top right`, `bottom left`, `bottom right`],
  }

  let {
    color_scale = $bindable(`interpolateViridis`),
    active_element = $bindable(null),
    models = $bindable(model_names_with_errors),
    current_model = $bindable([models[0] ?? ``]),
    manual_cbar_max = $bindable(false),
    normalized = $bindable(true),
    cbar_max = $bindable(0.3),
    ...rest
  }: ComponentProps<typeof PeriodicTable> & {
    color_scale?: D3InterpolateName | ((num: number) => string)
    active_element?: ChemicalElement | null
    models?: string[]
    // Must be string[] instead of string for svelte-multiselect to be correctly restored by snapshot
    current_model?: string[]
    manual_cbar_max?: boolean
    normalized?: boolean
    cbar_max?: number | null
  } = $props()

  const test_set_std = each_errors[`Test set standard deviation`]

  // selected models resolved in selection order (drives segment order in split tiles)
  let selected_models = $derived.by(() => {
    const resolved = current_model
      .map((name) => models_with_errors.find((model) => model.model_name === name))
      .filter((model) => model !== undefined)
    return resolved.length > 0 ? resolved : models_with_errors.slice(0, 1)
  })

  // null = no data (null error, or null std in normalized mode). Tiles coerce null to
  // 0 (rendered as missing_color); the hover inset shows n/a instead so missing data
  // can't read as a perfect score
  const norm_error = (model_key: string, element: string): number | null => {
    const val = each_errors[model_key]?.[element]
    const denom = normalized ? test_set_std[element] : 1
    return val == null || !denom ? null : val / denom
  }

  // single model: scalar per element; multiple: array per element -> split tiles
  let heatmap_values = $derived(
    Object.fromEntries(
      Object.keys(test_set_std).map((element) => {
        const errors = selected_models.map(
          (model) => norm_error(model.model_key, element) ?? 0,
        )
        return [element, errors.length === 1 ? errors[0] : errors]
      }),
    ) as Record<ElementSymbol, number | number[]>,
  )
  let current_data_max = $derived(Math.max(...Object.values(heatmap_values).flat()))
  let cs_range = $derived<[number, number]>([
    0,
    manual_cbar_max ? (cbar_max ?? 0) : current_data_max,
  ])
  // selected_models (not current_model) so the title tracks the fallback model when
  // the bound selection doesn't resolve, e.g. a stale snapshot-restored model name
  let cbar_title = $derived(
    selected_models.length === 1
      ? `${selected_models[0].model_name} (${normalized ? `normalized` : `eV/atom`})`
      : `Element-projected error (${normalized ? `normalized` : `eV/atom`})`,
  )

  const snapshot_data = () => ({
    color_scale,
    current_model,
    cbar_max,
    manual_cbar_max,
    normalized,
  })
  export const snapshot = {
    capture: snapshot_data,
    restore: (values: ReturnType<typeof snapshot_data>) =>
      ({ color_scale, current_model, cbar_max, manual_cbar_max, normalized } = values),
  }
</script>

<p>
  This periodic table heatmap shows the MAE of model-predicted convex hull distance
  projected onto each element. The errors for every structure in the test set are
  projected onto the fraction of each element in the composition and averaged over all
  structures. The error is the absolute difference per atom between predicted and actual
  energy distance to the convex hull. Select up to 4 models to compare them side by side:
  each element tile splits into one segment per model.
</p>

<Select bind:selected={current_model} options={models} maxSelect={4} minSelect={1} />

{#if selected_models.length > 1}
  <div class="split-legend">
    {#each selected_models as model, idx (model.model_key)}
      <span>
        <strong>{model.model_name}</strong>
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
  {color_scale}
  bind:active_element
  color_scale_range={cs_range}
  tile_props={{ float_fmt: `.2` }}
  show_photo={false}
  missing_color="rgba(255,255,255,0.3)"
  {...rest}
>
  {#snippet inset()}
    {#if active_element}
      <TableInset style="align-content: center">
        {#if selected_models.length === 1}
          <PtableInset
            element={active_element}
            elem_counts={heatmap_values as Record<ElementSymbol, number>}
            show_percent={false}
            unit="<small style='font-weight: lighter;'>eV / atom</small>"
          />
        {:else}
          <strong class="multi-model-inset">
            {active_element.name}: {#each selected_models as model, idx (model.model_key)}
              {#if idx > 0}&ensp;{/if}
              {@const elem_error = norm_error(model.model_key, active_element.symbol)}
              <span>
                {model.model_name}
                <b>{elem_error == null ? `n/a` : format_num(elem_error)}</b>
              </span>
            {/each}
          </strong>
        {/if}
        <ColorBar
          title={cbar_title}
          title_side="top"
          {color_scale}
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
  .multi-model-inset {
    display: block;
    text-align: center;
    min-height: 18pt;
  }
  .multi-model-inset span {
    font-weight: lighter;
  }
</style>
