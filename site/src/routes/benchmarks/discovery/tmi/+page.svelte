<script lang="ts">
  import elem_prev from '$figs/element-prevalence-vs-error.jsonl'
  import ModelSelect from '$lib/ModelSelect.svelte'
  import { wide_legend } from '$lib/fig-helpers'
  import { UrlModelSelection } from '$lib/model-selection.svelte'
  import { bind_url_params } from '$lib/url-state.svelte'
  import { ScatterPlot } from 'matterviz/plot'
  import DiscoveryMetricFigs from './discovery-metric-figs.md'
  import ElementErrorsPtableHeatmap from './ElementErrorsPtableHeatmap.svelte'

  const elem_prev_options = elem_prev.models.map(({ model_key, label }) => ({
    label,
    value: model_key,
  }))
  const elem_prev_selection = new UrlModelSelection(() => ({
    options: elem_prev_options,
    defaults: elem_prev_options.slice(0, 3).map(({ value }) => value),
  }))
  bind_url_params(
    (params) => {
      elem_prev_selection.read(params)
      if (!elem_prev_selection.selected.length) {
        elem_prev_selection.selected = elem_prev_options.slice(0, 3)
      }
    },
    () => [elem_prev_selection.url_entry],
  )

  const elem_prev_models = $derived(
    elem_prev.models.filter(({ model_key }) =>
      elem_prev_selection.values.includes(model_key),
    ),
  )
  const numeric_pairs = (
    x_values: (number | null)[],
    y_values: (number | null)[],
    elements: string[],
  ) => {
    const numeric_x_values: number[] = []
    const numeric_y_values: number[] = []
    const metadata: { elem: string }[] = []
    for (const [idx, x_val] of x_values.entries()) {
      const y_val = y_values[idx]
      if (x_val == null || y_val == null) continue
      numeric_x_values.push(x_val)
      numeric_y_values.push(y_val)
      metadata.push({ elem: elements[idx] })
    }
    return { x: numeric_x_values, y: numeric_y_values, metadata }
  }
</script>

<h1 id="discovery-too-much-information">Discovery: Too Much Information</h1>

Discovery diagnostics that didn't make the cut into the
<a href="/benchmarks/discovery">task page</a>.

<h2 id="per-element-model-error-heatmaps" style="text-align: center">
  Per-Element Model Error Heatmaps
</h2>

<ElementErrorsPtableHeatmap />

<br />

<DiscoveryMetricFigs />

<h2
  id="does-error-correlate-with-element-prevalence-in-training-set"
  style="text-align: center"
>
  Does error correlate with element prevalence in training set?
</h2>

Answer: not much. You might expect the more examples of structures containing a certain
element models have seen in the training set, the smaller their average error on test set
structures containing that element. That's not what we see in this plot. E<sub
  >above hull</sub
>
is all over the place as a function of elemental training set prevalence. Could be because the
error is dominated by the least abundant element in composition or the model errors are more
dependent on geometry than chemistry.

<label>
  Models
  <ModelSelect
    options={elem_prev_options}
    bind:value={elem_prev_selection.selected}
    min_select={1}
  />
</label>
<ScatterPlot
  series={elem_prev_models.map(({ label, color, y: error_values }) => ({
    ...numeric_pairs(elem_prev.occurrences, error_values, elem_prev.elements),
    label,
    markers: `points` as const,
    point_style: { fill: color },
  }))}
  x_axis={{ label: `MP Occurrences`, range: [0, null], format: `~s` }}
  y_axis={{ label: `Error (eV/atom)` }}
  legend={wide_legend}
>
  {#snippet tooltip({ x_formatted, y_formatted, metadata, label })}
    <strong>{metadata?.elem}</strong> ({label})<br />
    {x_formatted} MP occurrences<br />
    error: {y_formatted} eV/atom
  {/snippet}
</ScatterPlot>
