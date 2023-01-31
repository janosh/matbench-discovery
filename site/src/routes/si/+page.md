<script lang="ts">
  import RunTimePie from '$figs/2023-01-26-model-run-times-pie.svelte'
  import RocModels from '$figs/2023-01-30-roc-models.svelte'
</script>

## Model Run Times

{#if typeof document !== `undefined`}
<RunTimePie />

<figcaption>
@label:fig:model-run-times-pie
Creating this benchmark (excluding debugging runs) used a total of 3137 hours of compute time (mix of CPU and GPU, mostly CPU). Notably, the vast majority of that was used in the Bayesian optimization step of the BOWSR+MEGnet model.
</figcaption>
{/if}

{#if typeof document !== `undefined`}
<RocModels  />

<figcaption>
  @label:fig:roc-models
  Receiver operating characteristic (ROC) curves for each model. TPR/FPR is the true and false positive rates. Points are colored by their stability threshold. A material is classified as stable if E<sub>above hull</sub> is less than the stability threshold. Since all models predict E<sub>form</sub> (and M3GNet predicted energies are converted to formation energy before stability classification), they are insensitive to changes in the threshold.
</figcaption>
{/if}
