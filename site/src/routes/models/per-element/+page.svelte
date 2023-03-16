<script lang="ts">
  import { ElemCountInset } from '$lib'
  import type { ChemicalElement } from 'elementari'
  import { ColorScaleSelect, PeriodicTable, TableInset } from 'elementari'
  import { MultiSelect } from 'svelte-multiselect'
  import per_elem_errors from './per-element-model-each-errors.json'

  let color_scale = [`Oranges`]
  let active_element: ChemicalElement
  // $: active_counts = elem_counts[filter]
  let models = Object.keys(per_elem_errors)
  let current_model: string[] = [models[1]]
  $: heatmap_values = Object.values(per_elem_errors[current_model[0]])

  export const snapshot = {
    capture: () => ({ color_scale, current_model }),
    restore: (values) => ({ color_scale, current_model } = values),
  }
</script>

<h1>Too Much Information</h1>

Stuff that didn't make the cut into the
<a href="/models">model page</a>.

<h2>Per-Element Model Error Heatmaps</h2>

This periodic table is shaded by the average model error for each element. The errors for
every structure in the test set are projected onto the fraction of each element in the
composition and averaged. The error is the absolute difference between the predicted and
actual energy distance to the convex hull.

<MultiSelect bind:selected={current_model} options={models} maxSelect={1} minSelect={1} />

<ColorScaleSelect bind:selected={color_scale} />
<PeriodicTable
  {heatmap_values}
  color_scale={color_scale[0]}
  bind:active_element
  style="margin: 0 -2em;"
>
  <TableInset slot="inset" style="place-content: center;">
    <ElemCountInset element={active_element} elem_counts={heatmap_values} />
  </TableInset>
</PeriodicTable>
