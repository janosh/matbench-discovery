<script lang="ts">
  import { PtableInset } from '$lib'
  import type { ChemicalElement } from 'elementari'
  import { ColorBar, ColorScaleSelect, PeriodicTable, TableInset } from 'elementari'
  import { MultiSelect } from 'svelte-multiselect'
  import per_elem_errors from './per-element-each-errors.json'

  export let color_scale: string[] = [`Inferno`]
  export let active_element: ChemicalElement | null = null
  // $: active_counts = elem_counts[filter]
  export let models: string[] = Object.keys(per_elem_errors)
  export let current_model: string[] = [models[2]]
  export let manual_cbar_max: boolean = false
  export let normalized: boolean = true
  export let cbar_max: number | null = 0.3

  const test_set_std_key = Object.keys(per_elem_errors).find((key) =>
    key.includes(`Test set standard deviation`)
  ) as string
  const test_set_std = Object.values(per_elem_errors[test_set_std_key]) as number[]

  $: heatmap_values = (Object.values(per_elem_errors[current_model[0]]) as number[]).map(
    (val, idx) => {
      const denom = normalized ? test_set_std[idx] : 1
      if (denom) return val / denom
      return null
    }
  )
  $: current_data_max = Math.max(...heatmap_values)
  $: cs_range = [0, manual_cbar_max ? cbar_max : current_data_max]

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

This periodic table is shaded by the MAE for the model-predicted convex hull distance for
each element. The errors for every structure in the test set are projected onto the
fraction of each element in the composition and averaged over all structures. The error is
the absolute difference per atom between predicted and actual energy distance to the
convex hull.

<MultiSelect bind:selected={current_model} options={models} maxSelect={1} minSelect={1} />

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
    Divide errors by test set energies std. dev. over structures containing each element
    <input type="checkbox" bind:checked={normalized} />
  </label>
</form>

<PeriodicTable
  {heatmap_values}
  color_scale={color_scale[0]}
  bind:active_element
  color_scale_range={cs_range}
>
  <TableInset slot="inset" style="align-content: center;">
    <PtableInset
      element={active_element}
      elem_counts={heatmap_values}
      show_percent={false}
      unit="<small style='font-weight: lighter;'>eV / atom</small>"
    />
    <ColorBar
      text="{current_model[0]} ({normalized ? `normalized` : `eV/atom`})"
      text_side="top"
      color_scale={color_scale[0]}
      tick_labels={5}
      range={cs_range}
      style="width: 85%; margin: 0 2em;"
    />
  </TableInset>
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
</style>
