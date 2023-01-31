<script lang="ts">
  import RunTimePie from '$figs/2023-01-26-model-run-times-pie.svelte'
  import RocModels from '$figs/2023-01-30-roc-models.svelte'
</script>

# Supplementary Information

## Model Run Times

{#if typeof document !== `undefined`}
<RunTimePie style="margin: 1em;" />

<figcaption>
@label:fig:model-run-times-pie
Creating this benchmark (excluding debugging runs) used a total of 3137 hours of compute time (mix of CPU and GPU, mostly CPU). Notably, the vast majority of that was used in the Bayesian optimization step of the BOWSR+MEGnet model.
</figcaption>
{/if}

## ROC Curves

{#if typeof document !== `undefined`}
<RocModels  />

<figcaption>
  @label:fig:roc-models
  Receiver operating characteristic (ROC) curve for each model. TPR/FPR is the true and false positive rates. Points are colored by their stability threshold. A material is classified as stable if E<sub>above hull</sub> is less than the stability threshold. Since all models predict E<sub>form</sub> (and M3GNet predicted energies are converted to formation energy before stability classification), they are insensitive to changes in the threshold.
</figcaption>
{/if}

## Formation Energy MAE = Hull Distance MAE

To avoid potential confusion for people reading the code, we often calculate the formation energy MAE and report it as the MAE for the energy above the convex hull prediction. The former is more easily calculated but the two quantities are really the same. The formation energy of a material is the difference in energy between a material and its constituent elements in their standard states. The distance to the convex hull is defined as the difference between a material's formation energy and the minimum formation energy of all possible stable materials made from the same elements. Since the formation energy of a material is used to calculate the distance to the convex hull, the error of a formation energy prediction directly determines the error in the distance to the convex hull prediction.
