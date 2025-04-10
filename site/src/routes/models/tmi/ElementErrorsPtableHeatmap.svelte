<script lang="ts">
  import per_elem_each_errors from '$figs/per-element-each-errors.json'
  import { MODELS, PtableInset } from '$lib'
  import type { ChemicalElement } from 'elementari'
  import { ColorBar, ColorScaleSelect, PeriodicTable, TableInset } from 'elementari'
  import Select from 'svelte-multiselect'

  interface Props {
    color_scale?: string[]
    active_element?: ChemicalElement | null
    models?: string[]
    // must be string[] instead of string for svelte-multiselect to be correctly restored by snapshot
    current_model?: string[]
    manual_cbar_max?: boolean
    normalized?: boolean
    cbar_max?: number | null
  }

  let {
    color_scale = $bindable([`Viridis`]),
    active_element = $bindable(null),
    models = $bindable(MODELS.map((m) => m.model_name)),
    current_model = $bindable([models[2]]),
    manual_cbar_max = $bindable(false),
    normalized = $bindable(true),
    cbar_max = $bindable(0.3),
  }: Props = $props()

  const test_set_std_key = `Test set standard deviation`

  const test_set_std = per_elem_each_errors[test_set_std_key]

  let heatmap_values = $derived(
    Object.entries(per_elem_each_errors[current_model[0]]).map(([key, val]) => {
      const denom = normalized ? test_set_std[key] : 1
      if (denom) return val / denom
      return null
    }),
  )
  let current_data_max = $derived(Math.max(...heatmap_values))
  let cs_range = $derived([0, manual_cbar_max ? cbar_max : current_data_max])

  export const snapshot = {
    capture: () => ({
      color_scale,
      current_model,
      cbar_max,
      manual_cbar_max,
      normalized,
    }),
    restore: (values) =>
      ({ color_scale, current_model, cbar_max, manual_cbar_max, normalized } = values),
  }
</script>

<p style="max-width: 45em; margin: auto;">
  This periodic table heatmap shows the MAE of model-predicted convex hull distance
  projected onto each element. The errors for every structure in the test set are
  projected onto the fraction of each element in the composition and averaged over all
  structures. The error is the absolute difference per atom between predicted and actual
  energy distance to the convex hull.
</p>

<Select bind:selected={current_model} options={models} maxSelect={1} minSelect={1} />

<ColorScaleSelect bind:selected={color_scale} />

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
    Divide each element value by its std. dev. of target energies over all test structures
    containing a given element
    <input type="checkbox" bind:checked={normalized} />
  </label>
  <small>
    This is meant to correct for the fact that some elements are inherently more difficult
    to predict since some have a wider distribution of convex hull distances.
  </small>
</form>

<PeriodicTable
  {heatmap_values}
  color_scale={color_scale[0]}
  bind:active_element
  color_scale_range={cs_range}
  tile_props={{ precision: `0.2` }}
  show_photo={false}
>
  {#snippet inset()}
    <TableInset style="align-content: center;">
      <PtableInset
        element={active_element}
        elem_counts={heatmap_values}
        show_percent={false}
        unit="<small style='font-weight: lighter;'>eV / atom</small>"
      />
      <ColorBar
        text="{current_model[0]} ({normalized ? `normalized` : `eV/atom`})"
        label_side="top"
        color_scale={color_scale[0]}
        tick_labels={5}
        range={cs_range}
        style="width: 85%; margin: 0 2em;"
      />
    </TableInset>
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
</style>
