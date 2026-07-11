<script lang="ts">
  import elem_prev from '$figs/element-prevalence-vs-error.jsonl'
  import hist_largest from '$figs/hist-largest-each-errors-fp-diff.jsonl'
  import each_errors from '$figs/scatter-largest-each-errors-fp-diff.jsonl'
  import fp_diff from '$figs/scatter-largest-fp-diff-each-error.jsonl'
  import { ModelSelect } from '$lib'
  import { dashed, plotly_blue, plotly_red, wide_legend } from '$lib/fig-helpers'
  import { bind_url_params, valid_query_param } from '$lib/url-state.svelte'
  import type { UrlParamEntry } from '$lib/url-state.svelte'
  import DiscoveryMetricFigs from '$routes/models/tmi/discovery-metric-figs.md'
  import ElementErrorsPtableHeatmap from '$routes/models/tmi/ElementErrorsPtableHeatmap.svelte'
  import { BarPlot, BinnedScatterPlot, ScatterPlot } from 'matterviz/plot'

  // payload models arrive pre-styled (stable MODELS colors + leaderboard order) from the
  // json_payload plugin, so each dropdown below defaults to a top model
  const fp_diff_label = `|SSFP<sub>initial</sub> - SSFP<sub>final</sub>|`

  // per-figure model selection via dropdowns (faster than the old all-series-behind-a-
  // huge-legend figs). bind the model labels (svelte-multiselect string options; binding
  // the model objects directly would JSON-serialize their arrays into the DOM via the
  // hidden form-validation input) and look up the matching model object(s) for plotting.
  const find_model = <T extends { label: string }>(models: T[], label: string): T =>
    models.find((model) => model.label === label) ?? models[0]

  const elem_prev_labels = elem_prev.models.map((model) => model.label)
  const default_elem_prev = elem_prev_labels.slice(0, 3)
  // per-figure single-model dropdowns: one URL param each, defaulting to the top model
  const single_selects: Record<string, string[]> = {
    fp_model: fp_diff.models.map((model) => model.label),
    each_model: each_errors.models.map((model) => model.label),
    hist_model: hist_largest.models.map((model) => model.label),
  }

  let elem_prev_selected = $state([...default_elem_prev])
  let picked = $state(
    Object.fromEntries(
      Object.entries(single_selects).map(([key, model_labels]) => [key, model_labels[0]]),
    ),
  )

  // serialize the multi-select in canonical payload order so URL comparison against
  // the default is insensitive to the order models were clicked in
  const elem_prev_param = (selected: string[]): string =>
    elem_prev_labels.filter((label) => selected.includes(label)).join(`,`)

  const read_url_params = (params: URLSearchParams) => {
    const parsed = params
      .get(`models`)
      ?.split(`,`)
      .filter((label) => elem_prev_labels.includes(label))
    elem_prev_selected = parsed?.length ? parsed : [...default_elem_prev]
    for (const [key, model_labels] of Object.entries(single_selects)) {
      picked[key] = valid_query_param(params, key, model_labels[0], new Set(model_labels))
    }
  }
  bind_url_params(read_url_params, () => [
    [`models`, elem_prev_param(elem_prev_selected), elem_prev_param(default_elem_prev)],
    ...Object.entries(single_selects).map(
      ([key, model_labels]): UrlParamEntry => [key, picked[key], model_labels[0]],
    ),
  ])

  const elem_prev_models = $derived(
    elem_prev.models.filter((model) => elem_prev_selected.includes(model.label)),
  )
  const fp_diff_active = $derived(find_model(fp_diff.models, picked.fp_model))
  const each_errors_active = $derived(find_model(each_errors.models, picked.each_model))
  const hist_largest_active = $derived(find_model(hist_largest.models, picked.hist_model))

  // x extent of the shared fingerprint-diff values for the MAE ref line
  const fp_diff_extent = [Math.min(...fp_diff.fp_diff), Math.max(...fp_diff.fp_diff)]

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
    options={elem_prev_labels}
    bind:selected={elem_prev_selected}
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
    bind:value={picked.fp_model}
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
  density={{ color_scale: { type: `log`, scheme: `interpolateMagma` }, color_bar: null }}
  overlays={{
    ref_lines: [
      {
        x1: fp_diff_extent[0],
        y1: fp_diff_active.mae,
        x2: fp_diff_extent[1],
        y2: fp_diff_active.mae,
        ...dashed,
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
    bind:value={picked.each_model}
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
    bind:value={picked.hist_model}
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
/>
