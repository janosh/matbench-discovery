<script lang="ts">
  import elem_prev from '$figs/element-prevalence-vs-error.json.gz'
  import hist_largest from '$figs/hist-largest-each-errors-fp-diff.json.gz'
  import each_errors from '$figs/scatter-largest-each-errors-fp-diff.json.gz'
  import fp_diff from '$figs/scatter-largest-fp-diff-each-error.json.gz'
  import { dashed, plotly_blue, plotly_red, wide_legend } from '$lib/fig-helpers'
  import { BarPlot, BinnedScatterPlot, ScatterPlot } from 'matterviz/plot'
  import Select from 'svelte-multiselect'
  import DiscoveryMetricFigs from './discovery-metric-figs.md'
  import ElementErrorsPtableHeatmap from './ElementErrorsPtableHeatmap.svelte'

  const fp_diff_label = `|SSFP<sub>initial</sub> - SSFP<sub>final</sub>|`

  // per-figure model selection via dropdowns (faster than the old all-series-behind-a-
  // huge-legend figs). bind the model labels (svelte-multiselect string options; binding
  // the model objects directly would JSON-serialize their arrays into the DOM via the
  // hidden form-validation input) and look up the matching model object(s) for plotting.
  const find_model = <T extends { label: string }>(models: T[], label: string): T =>
    models.find((mdl) => mdl.label === label) ?? models[0]

  let elem_prev_selected = $state(elem_prev.models.slice(0, 3).map((mdl) => mdl.label))
  let fp_diff_model = $state(fp_diff.models[0].label)
  let each_errors_model = $state(each_errors.models[0].label)
  let hist_largest_model = $state(hist_largest.models[0].label)

  const elem_prev_models = $derived(
    elem_prev.models.filter((mdl) => elem_prev_selected.includes(mdl.label)),
  )
  const fp_diff_active = $derived(find_model(fp_diff.models, fp_diff_model))
  const each_errors_active = $derived(find_model(each_errors.models, each_errors_model))
  const hist_largest_active = $derived(
    find_model(hist_largest.models, hist_largest_model),
  )

  // x extent of the shared fingerprint-diff values for the MAE ref line
  const fp_diff_extent = [Math.min(...fp_diff.fp_diff), Math.max(...fp_diff.fp_diff)]
</script>

<h1>Too Much Information</h1>

Stuff that didn't make the cut into the&nbsp;<a href="/models">model page</a>.

<h2>Per-Element Model Error Heatmaps</h2>

<ElementErrorsPtableHeatmap />

<br />

<DiscoveryMetricFigs />

<h2>Does error correlate with element prevalence in training set?</h2>

Answer: not much. You might expect the more examples of structures containing a certain
element models have seen in the training set, the smaller their average error on test set
structures containing that element. That's not what we see in this plot. E<sub>above
  hull</sub>
is all over the place as a function of elemental training set prevalence. Could be because
the error is dominated by the least abundant element in composition or the model errors
are more dependent on geometry than chemistry.

<label class="model-select">
  Models
  <Select
    options={elem_prev.models.map((mdl) => mdl.label)}
    bind:selected={elem_prev_selected}
    minSelect={1}
  />
</label>
<ScatterPlot
  series={elem_prev_models.map(({ label, color, y }) => ({
    x: elem_prev.occurrences,
    y,
    label,
    markers: `points` as const,
    point_style: { fill: color },
    // one point per element -> element symbol in the tooltip
    metadata: elem_prev.elements.map((elem) => ({ elem })),
  }))}
  x_axis={{ label: `MP Occurrences`, range: [0, null], format: `~s` }}
  y_axis={{ label: `Error (eV/atom)` }}
  legend={wide_legend}
  style="height: 440px; margin: 2em 0"
>
  {#snippet tooltip({ x_formatted, y_formatted, metadata, label })}
    <strong>{metadata?.elem}</strong> ({label})<br />
    {x_formatted} MP occurrences<br />
    error: {y_formatted} eV/atom
  {/snippet}
</ScatterPlot>

<h2>Does error correlate with relaxation change?</h2>

Taking structures with the largest difference in atomic environments before vs after
relaxation as measured by<code>matminer</code>'s
<a
  href="https://hackingmaterials.lbl.gov/matminer/matminer.featurizers.structure.html#matminer.featurizers.structure.sites.SiteStatsFingerprint"
>
  <code>SiteStatsFingerprint</code>
</a>
(which is volume independent so changes in fingerprint require ion migration or similar)
and plotting against that the absolute E<sub>above hull</sub> errors for each model.

<label class="model-select">
  Model
  <Select
    options={fp_diff.models.map((mdl) => mdl.label)}
    bind:value={fp_diff_model}
    minSelect={1}
    maxSelect={1}
  />
  <small>MAE = {fp_diff_active.mae} eV/atom (dashed line)</small>
</label>
<BinnedScatterPlot
  series={[{
    x: fp_diff.fp_diff,
    y: fp_diff_active.y,
    label: fp_diff_active.label,
    color: fp_diff_active.color,
  }]}
  x_axis={{ label: fp_diff_label }}
  y_axis={{ label: `|E<sub>above hull</sub> error| (eV/atom)` }}
  density={{ color_scale: { type: `log`, scheme: `interpolateMagma` }, color_bar: null }}
  overlays={{
    ref_lines: [{
      x1: fp_diff_extent[0],
      y1: fp_diff_active.mae,
      x2: fp_diff_extent[1],
      y2: fp_diff_active.mae,
      ...dashed,
    }],
  }}
  style="height: 440px; margin: 2em 0"
/>

Same plot except taking the structures with largest difference in atomic environments
(again measured by
<code>SiteStatsFingerprint</code> before vs after relaxation) and plotting all model
errors.

<label class="model-select">
  Model
  <Select
    options={each_errors.models.map((mdl) => mdl.label)}
    bind:value={each_errors_model}
    minSelect={1}
    maxSelect={1}
  />
  <small>MAE = {each_errors_active.mae} eV/atom</small>
</label>
<ScatterPlot
  series={[{
    x: each_errors_active.x,
    y: each_errors_active.y,
    label: each_errors_active.label,
    markers: `points` as const,
  }]}
  x_axis={{ label: fp_diff_label, range: [0, null] }}
  y_axis={{ label: `Absolute error (eV/atom)` }}
  legend={null}
  style="height: 440px; margin: 2em 0"
/>

Another way to plot this is as a histogram. This shows the difference in
SiteStatsFingerprint before vs after relaxation for structures with the largest (err<sub
>max</sub>) and smallest (err<sub>min</sub>) absolute error in predicted E<sub>above
  hull</sub>
for each model and the mean of all models.

<label class="model-select">
  Model
  <Select
    options={hist_largest.models.map((mdl) => mdl.label)}
    bind:value={hist_largest_model}
    minSelect={1}
    maxSelect={1}
  />
</label>
<BarPlot
  series={[
    { ...hist_largest_active.err_min, label: `err<sub>min</sub>`, color: plotly_blue },
    { ...hist_largest_active.err_max, label: `err<sub>max</sub>`, color: plotly_red },
  ]}
  mode="overlay"
  x_axis={{ label: fp_diff_label, range: [0, null] }}
  y_axis={{ label: `Count` }}
  show_legend
  show_controls={false}
  style="height: 440px; margin: 2em 0"
/>

<style>
  .model-select {
    display: flex;
    gap: 1ex;
    align-items: center;
    margin-top: 1em;
  }
  .model-select :global(.multiselect) {
    min-width: 16em;
  }
</style>
