<script lang="ts">
  import { onMount } from 'svelte'
  import BoxHullDistErrors from '$figs/box-hull-dist-errors.svelte'
  import CumulativePrecisionRecall from '$figs/cumulative-precision-recall.svelte'
  import EachParityModels from '$figs/each-parity-models-6x3.svelte'
  import HistClfPredHullDistModels from '$figs/hist-clf-pred-hull-dist-models-6x3.svelte'
  import RocModels from '$figs/roc-models.svelte'
  import RollingMaeVsHullDistModels from '$figs/rolling-mae-vs-hull-dist-models.svelte'

  let mounted: boolean = false
  onMount(() => (mounted = true))
</script>

{#if mounted}
<CumulativePrecisionRecall />

> @label:fig:cumulative-precision-recall Model precision and recall for thermodynamic stability classification as a function of number of materials ranked from most to least stable by each model.
> CHGNet initially achieves the highest cumulative precision and recall.
> Simulates materials discovery efforts of different sizes since a typical campaign will rank hypothetical materials by model-predicted hull distance from most to least stable and validate the most stable predictions first.
> A higher fraction of correct stable predictions corresponds to higher precision and fewer stable materials overlooked correspond to higher recall.
> This figure highlights how different models perform better or worse depending on the length of the discovery campaign.
> The UIPs (CHGNet, M3GNet, MACE) are seen to offer significantly improved precision on shorter campaigns of ~20k or less materials validated as they are less prone to false positive predictions among highly stable materials.

{/if}

{#if mounted}
<RocModels />

> @label:fig:roc-models Receiver operating characteristic (ROC) curve for each model. TPR/FPR = true/false positive rate. FPR on the $x$ axis is the fraction of unstable structures classified as stable. TPR on the $y$ axis is the fraction of stable structures classified as stable.

{/if}

{#if mounted}
<BoxHullDistErrors />

> @label:fig:box-hull-dist-errors Box plot of interquartile ranges (IQR) of hull distance errors for each model. The whiskers extend to the 5th and 95th percentiles. The horizontal line inside the box shows the median. BOWSR has the highest median error, while Voronoi RF has the highest IQR. Note that MEGNet and CGCNN are the only models with a positive median. Their hull distance errors are biased towards more frequently predicting thermodynamic instability, explaining why they are closest to getting the overall number of stable structures in the test set right (see cumulative precision/recall in @fig:cumulative-precision-recall).

{/if}

{#if mounted}
<RollingMaeVsHullDistModels style="place-self: center;" />

> @label:fig:rolling-mae-vs-hull-dist-models Universal potentials are more reliable classifiers because they exit the red triangle earliest.
> These lines show the rolling MAE on the WBM test set as the energy to the convex hull of the MP training set is varied.
> Lower is better.
> Inside the large red 'triangle of peril', models are most likely to misclassify structures.
> As long as a model's rolling MAE remains inside the triangle, its mean error is larger than the distance to the convex hull.
> If the model's error for a given prediction happens to point towards the stability threshold at $E$<sub>above MP hull</sub> = 0, its average error will change the stability classification from true positive/negative to false negative/positive.
> The width of the 'rolling window' box indicates the width over which prediction errors were averaged.

{/if}

{#if mounted}
<EachParityModels />

> @label:fig:each-parity-models Parity plots of model-predicted energy distance to the convex hull (based on their formation energy predictions) vs DFT ground truth, color-coded by log density of points.
> Models are sorted left to right and top to bottom by MAE.

{/if}

{#if mounted}
<HistClfPredHullDistModels />

> @label:fig:hist-clf-pred-hull-dist-models Distribution of model-predicted hull distance colored by stability classification. Models are sorted from top to bottom by F1 score. The thickness of the red and yellow bands shows how often models misclassify as a function of how far away from the convex hull they place a material. While CHGNet and M3GNet perform almost equally well overall, these plots reveal that they do so via different trade-offs. M3GNet commits fewer false negatives but more false positives predictions compared to CHGNet. In a real discovery campaign, false positives have a higher opportunity cost than false negatives, since they result in wasted DFT relaxations or even synthesis time in the lab. A false negative by contrast is just one missed opportunity out of many. For this reason, models with high true positive rate (TPR) even at the expense of lower true negative rate (TNR) are generally preferred.

{/if}
