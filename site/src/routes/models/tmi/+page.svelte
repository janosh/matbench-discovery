<script lang="ts">
  import { browser } from '$app/environment'
  import LargestEachErrorsFpDiffModels from '$figs/largest-each-errors-fp-diff-models.svelte'
  import LargeFpDiffVsEachError from '$figs/largest-fp-diff-each-error-models.svelte'
  import { PtableInset } from '$lib'
  import type { ChemicalElement } from 'elementari'
  import { ColorBar, ColorScaleSelect, PeriodicTable, TableInset } from 'elementari'
  import { MultiSelect } from 'svelte-multiselect'
  import per_elem_errors from './per-element-model-each-errors.json'

  let color_scale = [`Inferno`]
  let active_element: ChemicalElement
  // $: active_counts = elem_counts[filter]
  let models = Object.keys(per_elem_errors)
  let current_model: string[] = [models[1]]
  let manual_cbar_max: boolean = false
  $: heatmap_values = Object.values(per_elem_errors[current_model[0]]) as number[]
  $: current_data_max = Math.max(...heatmap_values)
  let cbar_max: number | null = 0.03

  export const snapshot = {
    capture: () => ({ color_scale, current_model, cbar_max, manual_cbar_max }),
    restore: (values) =>
      ({ color_scale, current_model, cbar_max, manual_cbar_max } = values),
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

<form>
  Manual color bar max
  <input type="checkbox" bind:checked={manual_cbar_max} />
  <input
    type="range"
    disabled={!manual_cbar_max}
    bind:value={cbar_max}
    min={0.01}
    max={0.15}
    step={0.001}
  />
  {cbar_max}
</form>

<PeriodicTable
  {heatmap_values}
  color_scale={color_scale[0]}
  bind:active_element
  style="margin: 0 -2em;"
  color_scale_range={[0, manual_cbar_max ? cbar_max : current_data_max]}
  precision={4}
>
  <TableInset slot="inset" style="align-content: center;">
    <PtableInset
      element={active_element}
      elem_counts={heatmap_values}
      precision={5}
      show_percent={false}
      style="min-height: 2em;"
      unit="<small style='font-weight: lighter;'>eV / atom</small>"
    />
    <ColorBar
      color_scale={color_scale[0]}
      tick_labels={5}
      range={[0, manual_cbar_max ? cbar_max : current_data_max]}
      precision={3}
      style="width: 100%; margin: 0 2em;"
    />
  </TableInset>
</PeriodicTable>

<h2>Analysis of E<sub>above hull</sub> Errors as a Function of Relaxation Change</h2>

Taking structures with the largest difference in atomic environments before vs after
relaxation as measured by<code>matminer</code>'s
<a
  href="https://hackingmaterials.lbl.gov/matminer/matminer.featurizers.structure.html#matminer.featurizers.structure.sites.SiteStatsFingerprint"
>
  <code>SiteStatsFingerprint</code>
</a>
(which is volume independent so changes in fingerprint require ion migration or similar) and
plotting against that the absolute E<sub>above hull</sub> errors for each model.

{#if browser}
  <LargeFpDiffVsEachError style="margin: 2em 0;" />
{/if}

Same plot except taking the structures with largest difference in atomic environments
(again measured by
<code>SiteStatsFingerprint</code> before vs after relaxation) and plotting all model
errors.

{#if browser}
  <LargestEachErrorsFpDiffModels style="margin: 2em 0;" />
{/if}

<style>
  form {
    display: flex;
    place-content: center;
    margin: 1em;
    gap: 1ex;
  }
</style>
