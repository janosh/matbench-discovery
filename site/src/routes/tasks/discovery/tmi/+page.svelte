<script lang="ts">
  import elem_prev from '$figs/element-prevalence-vs-error.jsonl'
  import hist_largest from '$figs/hist-largest-each-errors-fp-diff.jsonl'
  import each_errors from '$figs/scatter-largest-each-errors-fp-diff.jsonl'
  import fp_diff from '$figs/scatter-largest-fp-diff-each-error.jsonl'
  import { ModelSelect } from '$lib'
  import { dashed, series_blue, series_red, wide_legend } from '$lib/fig-helpers'
  import { UrlModelSelection } from '$lib/model-selection.svelte'
  import { bind_url_params, valid_query_param } from '$lib/url-state.svelte'
  import type { UrlParamEntry } from '$lib/url-state.svelte'
  import { BarPlot, BinnedScatterPlot, ScatterPlot } from 'matterviz/plot'
  import DiscoveryMetricFigs from './discovery-metric-figs.md'
  import ElementErrorsPtableHeatmap from './ElementErrorsPtableHeatmap.svelte'

  // payload models arrive pre-styled (stable MODELS colors + leaderboard order) from the
  // json_payload plugin, so each dropdown below defaults to a top model
  const fp_diff_label = `|SSFP<sub>initial</sub> - SSFP<sub>final</sub>|`

  // Selections and URL state use immutable model keys; labels remain display-only.
  const find_model = <T extends { model_key: string }>(
    models: T[],
    model_key: string,
  ): T => models.find((model) => model.model_key === model_key) ?? models[0]
  type ModelOption = { label: string; value: string }
  const model_options = (models: { model_key: string; label: string }[]): ModelOption[] =>
    models.map(({ model_key, label }) => ({ label, value: model_key }))

  const elem_prev_options = model_options(elem_prev.models)
  const elem_prev_selection = new UrlModelSelection(() => ({
    options: elem_prev_options,
    defaults: elem_prev_options.slice(0, 3).map(({ value }) => value),
  }))
  // per-figure single-model dropdowns: one URL param each, defaulting to the top model
  const single_selects = {
    fp_model: model_options(fp_diff.models),
    each_model: model_options(each_errors.models),
    hist_model: model_options(hist_largest.models),
  }

  let picked = $state<Record<string, ModelOption[]>>(
    Object.fromEntries(
      Object.entries(single_selects).map(([key, options]) => [key, [options[0]]]),
    ),
  )

  const read_url_params = (params: URLSearchParams) => {
    elem_prev_selection.read(params)
    if (!elem_prev_selection.selected.length) {
      elem_prev_selection.selected = elem_prev_options.slice(0, 3)
    }
    for (const [key, options] of Object.entries(single_selects)) {
      const options_by_key = Object.fromEntries(
        options.map((option) => [option.value, option]),
      )
      const model_key = valid_query_param(params, key, options[0].value, options_by_key)
      picked[key] = [options_by_key[model_key]]
    }
  }
  bind_url_params(read_url_params, () => [
    elem_prev_selection.url_entry,
    ...Object.entries(single_selects).map(([key, options]): UrlParamEntry => [
      key,
      picked[key][0].value,
      options[0].value,
    ]),
  ])

  const elem_prev_models = $derived(
    elem_prev.models.filter(({ model_key }) =>
      elem_prev_selection.values.includes(model_key),
    ),
  )
  const fp_diff_active = $derived(find_model(fp_diff.models, picked.fp_model[0].value))
  const each_errors_active = $derived(
    find_model(each_errors.models, picked.each_model[0].value),
  )
  const hist_largest_active = $derived(
    find_model(hist_largest.models, picked.hist_model[0].value),
  )

  const numeric_pairs = (
    x_values: (number | null)[],
    y_values: (number | null)[],
    elements: string[] = [],
  ) => {
    const numeric_x_values: number[] = []
    const numeric_y_values: number[] = []
    const metadata: { elem?: string }[] = []
    for (const [idx, x_val] of x_values.entries()) {
      const y_val = y_values[idx]
      if (x_val == null || y_val == null) continue
      numeric_x_values.push(x_val)
      numeric_y_values.push(y_val)
      if (elements.length) metadata.push({ elem: elements[idx] })
    }
    const numeric_values = { x: numeric_x_values, y: numeric_y_values }
    return elements.length ? { ...numeric_values, metadata } : numeric_values
  }
</script>

<h1>Discovery: Too Much Information</h1>

Discovery diagnostics that didn't make the cut into the
<a href="/tasks/discovery">task page</a>.

<h2>Per-Element Model Error Heatmaps</h2>

<ElementErrorsPtableHeatmap />

<br />

<DiscoveryMetricFigs />

<h2>Does error correlate with element prevalence in training set?</h2>

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
    bind:selected={elem_prev_selection.selected}
    minSelect={1}
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

<h2>Does error correlate with relaxation change?</h2>

Taking structures with the largest difference in atomic environments before vs after
relaxation as measured by<code>matminer</code>'s
<a
  href="https://hackingmaterials.lbl.gov/matminer/matminer.featurizers.structure.html#matminer.featurizers.structure.sites.SiteStatsFingerprint"
>
  <code>SiteStatsFingerprint</code>
</a>
(which is volume independent so changes in fingerprint require ion migration or similar) and
plotting against that the absolute E<sub>above hull</sub> errors for each model.

<label>
  Model
  <ModelSelect
    options={single_selects.fp_model}
    bind:selected={picked.fp_model}
    minSelect={1}
    maxSelect={1}
  />
  <small>MAE = {fp_diff_active.mae} eV/atom (dashed line)</small>
</label>
<BinnedScatterPlot
  series={[
    {
      ...numeric_pairs(fp_diff.fp_diff, fp_diff_active.y),
      label: fp_diff_active.label,
      color: fp_diff_active.color,
    },
  ]}
  x_axis={{ label: fp_diff_label }}
  y_axis={{ label: `|E<sub>above hull</sub> error| (eV/atom)` }}
  density={{ color_scale: { type: `log`, scheme: `interpolateMagma` } }}
  color_bar={null}
  overlays={{
    ref_lines: [
      {
        type: `horizontal`,
        y: fp_diff_active.mae,
        style: dashed,
      },
    ],
  }}
/>

Same plot except taking the structures with largest difference in atomic environments
(again measured by
<code>SiteStatsFingerprint</code> before vs after relaxation) and plotting all model
errors.

<label>
  Model
  <ModelSelect
    options={single_selects.each_model}
    bind:selected={picked.each_model}
    minSelect={1}
    maxSelect={1}
  />
  <small>MAE = {each_errors_active.mae} eV/atom</small>
</label>
<ScatterPlot
  series={[
    {
      x: each_errors_active.x,
      y: each_errors_active.y,
      label: each_errors_active.label,
      markers: `points` as const,
    },
  ]}
  x_axis={{ label: fp_diff_label, range: [0, null] }}
  y_axis={{ label: `Absolute error (eV/atom)` }}
  legend={null}
/>

Another way to plot this is as a histogram. This shows the difference in
SiteStatsFingerprint before vs after relaxation for structures with the largest (err<sub
  >max</sub
>) and smallest (err<sub>min</sub>) absolute error in predicted E<sub>above hull</sub> for
each model and the mean of all models.

<label>
  Model
  <ModelSelect
    options={single_selects.hist_model}
    bind:selected={picked.hist_model}
    minSelect={1}
    maxSelect={1}
  />
</label>
<BarPlot
  series={[
    { ...hist_largest_active.err_min, label: `err<sub>min</sub>`, color: series_blue },
    { ...hist_largest_active.err_max, label: `err<sub>max</sub>`, color: series_red },
  ]}
  mode="overlay"
  x_axis={{ label: fp_diff_label, range: [0, null] }}
  y_axis={{ label: `Count` }}
  show_legend
  show_controls={false}
/>
