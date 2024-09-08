<script lang="ts">
  import { browser } from '$app/environment'
  import EachErrorVsLeastPrevalentElementInStruct from '$figs/each-error-vs-least-prevalent-element-in-struct.svelte'
  import ElementPrevalenceVsErr from '$figs/element-prevalence-vs-error.svelte'
  import HistLargestEachErrorsFpDiffModels from '$figs/hist-largest-each-errors-fp-diff-models.svelte'
  import ScatterLargestEachErrorsFpDiffModels from '$figs/scatter-largest-each-errors-fp-diff-models.svelte'
  import ScatterLargeFpDiffVsEachError from '$figs/scatter-largest-fp-diff-each-error-models.svelte'
  import ElementErrorsPtableHeatmap from './ElementErrorsPtableHeatmap.svelte'
</script>

<h1>Too Much Information</h1>

Stuff that didn't make the cut into the&nbsp;<a href="/models">model page</a>.

<h2 style="margin: 4em auto 1em; text-align: center;">
  Per-Element Model Error Heatmaps
</h2>

<ElementErrorsPtableHeatmap />

<h2>Does error correlate with element prevalence in training set?</h2>

Answer: not much. You might (or might not) expect the more examples of structures
containing a certain element models have seen in the training set, the smaller their
average error on test set structures containing that element. That's not what we see in
this plot. E<sub>above hull</sub> is all over the place as a function of elemental
training set prevalence. Could be because the error is dominated by the least abundant
element in composition or the model errors are more dependent on geometry than chemistry.

{#if browser}
  <ElementPrevalenceVsErr style="margin: 2em 0;" />
{/if}

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

{#if browser}
  <ScatterLargeFpDiffVsEachError style="margin: 2em 0;" />
{/if}

Same plot except taking the structures with largest difference in atomic environments
(again measured by
<code>SiteStatsFingerprint</code> before vs after relaxation) and plotting all model
errors.

{#if browser}
  <ScatterLargestEachErrorsFpDiffModels style="margin: 2em 0;" />
{/if}

Another way to plot this is as a histogram. This shows the difference in
SiteStatsFingerprint before vs after relaxation for structures with the largest (err<sub
  >max</sub
>) and smallest (err<sub>min</sub>) absolute error in predicted E<sub>above hull</sub>
for each model and the mean of all models.

{#if browser}
  <HistLargestEachErrorsFpDiffModels style="margin: 2em 0;" />
{/if}

<h2>
  Does model error correlate with structure's least prevalent element in training set?
</h2>

Answer: a little. The fact that structures containing only elements with high prevalence
in the training set consistently see low errors across all models suggests that to get a
universally robust model, it needs to be trained on lots of examples for every element.

{#if browser}
  <EachErrorVsLeastPrevalentElementInStruct style="margin: 2em 0;" />
{/if}
