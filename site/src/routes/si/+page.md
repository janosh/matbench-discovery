<script lang="ts">
  import MetricsTableMegnetUipCombos from '$figs/metrics-table-uip-megnet-combos.svelte'
  import MetricsTableFirst10k from '$figs/metrics-table-first-10k.svelte'
  import RunTimeBars from '$figs/model-run-times-bar.svelte'
  import CumulativeMaeRmse from '$figs/cumulative-mae-rmse.svelte'
  import BoxHullDistErrors from '$figs/box-hull-dist-errors.svelte'
  import RocModels from '$figs/roc-models-all-in-one.svelte'
  import { browser } from '$app/environment'
  import MPRefEnergies from '$figs/mp-elemental-ref-energies.svelte'
  import WrenformerRollingMaeBatches from '$figs/rolling-mae-vs-hull-dist-wbm-batches-wrenformer.svelte'
  import CHGNetRollingMaeBatches from '$figs/rolling-mae-vs-hull-dist-wbm-batches-chgnet.svelte'
  import M3gnetRollingMaeBatches from '$figs/rolling-mae-vs-hull-dist-wbm-batches-m3gnet.svelte'
  import MegnetRollingMaeBatches from '$figs/rolling-mae-vs-hull-dist-wbm-batches-megnet.svelte'
  import CgcnnRollingMaeBatches from '$figs/rolling-mae-vs-hull-dist-wbm-batches-cgcnn.svelte'
  import VoronoiRfRollingMaeBatches from '$figs/rolling-mae-vs-hull-dist-wbm-batches-voronoi-rf.svelte'
  import SpacegroupSunburstWbm from '$figs/spacegroup-sunburst-wbm.svelte'
  import SpacegroupSunburstWrenformerFailures from '$figs/spacegroup-sunburst-wrenformer-failures.svelte'
  import LargestErrorScatterSelect from './largest-error-scatter-select.svelte'
  import HullDistScatterWrenformerFailures from '$figs/hull-dist-scatter-wrenformer-failures.svelte'
  import ProtoCountsWrenformerFailures from '$figs/proto-counts-wrenformer-failures.svelte'
  import ElementPrevalenceVsError from '$figs/element-prevalence-vs-error.svelte'
  import { onMount } from 'svelte'

  let mounted = false
  onMount(() => (mounted = true))
</script>

# Supplementary Information

## Metrics for 10k materials predicted most stable

<MetricsTableFirst10k />

> @label:fig:metrics-table-first-10k An actual discovery campaign is unlikely to validate all stable predictions coming from a given model like we did in the [metrics table](/preprint#fig:metrics-table). Presumably, it will rank model predictions from most to least stable and go down that list as far time and compute budgets permit. Assuming that increases in compute resources will allow average future discovery campaigns to grow in scope, we believe 10 k model validations to be a reasonable cutoff. To simulate this scenario, we calculated classification and regression metrics for the 10 k test set materials predicted to be most stable by each model.<br>
> We again show dummy performance in the bottom row. Note that each model is now evaluated on a different slice of the data, but this is still dummy performance across the whole dataset. CHGNet and M3GNet achieve a very impressive 83% and 80% precision, respectively. In concrete terms, this means in a discovery campaign that validates 10 k model predictions from a search pool of 257 k crystals that are chemically dissimilar from the training set and of which 16.7 % are stable, CHGNet and M3GNet would deliver 4 stable structures for every 5 predictions validated. In light of the significant resulting increase in stability hit rate, these models are well worth integrating into future materials searches.

## ROC Curves

{#if mounted}
<RocModels />
{/if}

> @label:fig:roc-models Receiver operating characteristic (ROC) curve for each model. TPR/FPR = true/false positive rate. FPR on the $x$-axis is the fraction of unstable structures classified as stable. TPR on the $y$-axis is the fraction of stable structures classified as stable. Points are colored by stability threshold $t$ which sweeps from $-0.4 \ \frac{\text{eV}}{\text{atom}} \leq t \leq 0.4 \ \frac{\text{eV}}{\text{atom}}$ above the hull. A material is classified as stable if the predicted E<sub>above hull</sub> lies below the stability threshold. Since all models predict E<sub>form</sub> (and M3GNet predicted energies are converted to formation energy before stability classification), they are insensitive to changes in the threshold $t$. M3GNet wins in area under curve (AUC) with 0.87, coming in 34% higher than the worst model Voronoi Random Forest. The diagonal 'No skill' line shows performance of a dummy model that randomly ranks material stability.

## Cumulative MAE + RMSE

{#if mounted}
<CumulativeMaeRmse />
{/if}

> @label:fig:cumulative-mae-rmse Cumulative mean absolute error (MAE) and root mean square error (RMSE) during a simulated discovery campaign. This figure expands on the [precision-recall figure](/preprint#fig:cumulative-precision-recall). The $x$-axis again shows number of materials sorted by model-predicted stability or 'campaign length'. This allows the reader to choose a cutoff point given their discovery campaign's resource constraints for validating model predictions and then read off the optimal model given those constraints.
> CHGNet achieves the lowest regression error profile, with a larger gap to the runner-up model M3GNet than in the precision-recall plots. This is likely due to the difference in TPR/TNR trade off between CHGNet and M3GNet. M3GNet has TNR = 0.80 vs CHGNet's TNR = 0.87. Higher TNR means lower FPR. Lower false positive rate means lower cumulative MAE and RMSE. Lines end when models stop predicting materials as stable, so these cumulative plots only contain model-predicted positive (stable) materials. Besides the high opportunity cost of false positives, this highlights another reason to prioritize low FPR in discovery models: lower error on the predictions of highest relevance.

## Model Run Times

{#if mounted}
<RunTimeBars style="margin: 1em;" />
{/if}

> @label:fig:model-run-times-pie Creating this benchmark (excluding debugging runs) used a total of 3210 hours of compute time (mix of CPU and GPU, mostly CPU). Notably, the vast majority of that (2705 h) was used in the Bayesian optimization step of BOWSR.<br>
> Some bars have two sections. The bottom shows training time, the upper test time. If there's only one section the model was not re-trained for this benchmark and the bar shows only the test time.

## Formation Energy MAE = Hull Distance MAE

To avoid potential confusion for people reading the code, we may in places calculate the formation energy MAE and report it as the MAE for the energy above the convex hull prediction. The former is more easily calculated but the two quantities are the same. The formation energy of a material is the difference in energy between a material and its constituent elements in their standard states. The distance to the convex hull is defined as the difference between a material's formation energy and the minimum formation energy of all possible stable materials made from the same elements. Since the formation energy of a material is used to calculate the distance to the convex hull, the error of a formation energy prediction directly determines the error in the distance to the convex hull prediction.

A further point of clarification: whenever we say distance to the convex hull we mean the signed distance that is positive for thermodynamically unstable materials above and negative for stable materials below the hull.

## MP Elemental Reference Energies

{#if mounted}
<MPRefEnergies />
{/if}

> @label:fig:mp-elemental-reference-energies WBM formation energies were calculated w.r.t. these Materials Project elemental reference energies ([queried on 2023-02-07](https://github.com/janosh/matbench-discovery/blob/-/data/mp/2023-02-07-mp-elemental-reference-entries.json.gz)). Marker size indicates the number of atoms in the reference structure. Hover points for details.

## WBM Batch Robustness as a Measure of Extrapolation Prowess

The figures below show the rolling MAE as a function of distance to the convex hull for each of the 5 WBM batches of elemental substitution. As a reminder, the WBM test set was generated in 5 successive batches, each applying another element replacement according to ISCD-mined chemical similarity scores to MP source structures and new stable crystals generated in the previous round. Naively, you'd expect model performance to degrade with increasing batch count as repeated substitutions should on average 'diffuse' deeper into uncharted regions of material space, requiring the model to extrapolate more. We observe this effect for some models but not others.

{#if mounted}

<div style="display: grid; grid-template-columns: 1fr 1fr; margin: 0 -1em 0 -4em;">
  <M3gnetRollingMaeBatches style="aspect-ratio: 1.2;" />
  <CHGNetRollingMaeBatches style="aspect-ratio: 1.2;" />
  <WrenformerRollingMaeBatches style="aspect-ratio: 1.2;" />
  <MegnetRollingMaeBatches style="aspect-ratio: 1.2;" />
  <VoronoiRfRollingMaeBatches style="aspect-ratio: 1.2;" />
  <CgcnnRollingMaeBatches style="aspect-ratio: 1.2;" />
</div>
{/if}

> @label:fig:rolling-mae-vs-hull-dist-wbm-batches-models On WBM batch 1, Wrenformer performs comparably to the SOTA UIPs M3GNet and CHGNet. However, the M3GNet and CHGNet UIPs show stronger extrapolative performance as they barely deteriorate in performance on later batches that move further away from the original training distribution. Wrenformer by contrast exhibits a pronounced increase in MAE with batch count.

MEGNet and CGCNN both incur a higher rolling MAE than Wrenformer across all 5 batches and across most or all of the hull distance range visible in these plots. However, similar to the UIPs, MEGNet and CGCNN exhibit very little degradation for higher batch counts. The fact that higher errors for later batches is specific to Wrenformer and Voronoi RF suggests that training only on composition and coarse-grained structural features (spacegroup and Wyckoff positions in the case of Wrenformer; coordination numbers, local environment properties, etc. in the case of Voronoi RF) alone is insufficient to learn an extrapolatable map of the PES.

Given its strong performance on batch 1, it is possible that given sufficiently diverse training data, Wrenformer could become similarly accurate to the UIPs across the whole PES landscape at substantially less training and inference cost. However, the loss of predictive forces and stresses may make Wrenformer unattractive for certain applications even then.

## Largest Errors vs DFT Hull Distance

{#if mounted}
<LargestErrorScatterSelect />
{/if}

> @label:fig:scatter-largest-errors-models-mean-vs-true-hull-dist DFT vs predicted hull distance (average over all models) for the 200 largest error structures colored by model disagreement (as measured by standard deviation in hull distance predictions from different models) and markers sized by number of atoms in the structures. This plot shows that high-error predictions are biased towards predicting too small hull distance. This is unsurprising considering MP training data mainly consists of low-energy structure.<br>
> However, note the clear color separation between the mostly blue low-energy-bias predictions and the yellow/red high error prediction. Blue means models are in good agreement, i.e. all models are "wrong" together. Red/yellow are large-error predictions with little model agreement, i.e. all models are wrong in different ways. It is possible that some of the blue points with large error yet good agreement among models are in fact accurate ML predictions for a DFT relaxation gone wrong. Zooming in on the blue points reveals that many of them are large. Larger markers correspond to larger structures where DFT failures are less surprising. This suggests ML model committees could be used to cheaply screen large databases for DFT errors in a high-throughput manner.

## MEGNet formation energies from UIP-relaxed structures

{#if mounted}
<MetricsTableMegnetUipCombos select={[`model`, `MEGNet`, `CHGNet`, `M3GNet`, `CHGNet + MEGNet`, `M3GNet + MEGNet`]} />
{/if}

> @label:fig:metrics-table-megnet-uip-combos This table shows metrics obtained by combining MEGNet with both UIPs. The metrics in rows labeled M3GNet→MEGNet and CHGNet→MEGNet are the result of passing M3GNet/CHGNet-relaxed structures into MEGNet for formation energy prediction. Both combos perform worse than using the respective UIPs on their own with a more pronounced performance drop from CHGNet to CHGNet→MEGNet than M3GNet to M3GNet→MEGnet. This suggests MEGNet has learned no additional knowledge of the PES that is not already present in the UIPs. However, both combos perform better than MEGNet on its own, demonstrating that UIP relaxation provides real utility at very low cost for any downstream structure-dependent analysis.

The UIPs M3GNet and CHGNet are both trained to predict DFT energies (including/excluding MP2020 energy corrections for CHGNet/M3GNet) while MEGNet is trained to predict formation energies.

The two types of target energies should in principle be equivalent since the difference between raw DFT energies and (corrected) formation energies is just a linear transformation of subtracting the elemental reference energies weighted by composition. In that sense, changing the target from DFT energies to formation energies merely constitutes a change of gauge. But in practice, one might expect the problem of stability prediction to be easier when trained on formation energies as the model output is one less step removed from the ultimate target of crystal stability. Empirical evidence so far suggests the opposite, however. The UIP + MEGNet combos both perform worse than UIPs by themselves, suggesting using formation energy as targets offers no benefits.

We highlight this here to refute the suggestion that training on raw DFT energies results in poorly calibrated predictions for deriving formation energies and subsequent stability predictions.

## Spacegroup prevalence in Wrenformer failure cases

{#if mounted}

<div style="display: flex; gap: 1em; justify-content: space-around; flex-wrap: wrap;">
<SpacegroupSunburstWrenformerFailures />
<SpacegroupSunburstWbm />
</div>
<HullDistScatterWrenformerFailures />
{/if}

> @label:fig:spacegroup-prevalence-wrenformer-failures The left spacegroup sunburst shows spacegroup 71 is by far the dominant lattice symmetry among the 941 Wrenformer failure cases where $E_\text{above hull,DFT} < 1$ and $E_\text{above hull,Wrenformer} > 1$ (points inside the shaded rectangle). On the right side for comparison is the spacegroup sunburst for the entire WBM test set.

Looking at the occurrence counts of isopointal prototypes in the shaded rectangle and comparing them with the occurrence of those same prototypes in the MP training data counts, we find almost no support for failing structure prototypes. This suggests the reason Wrenformer fails so spectacularly on these structures is that it cannot deal with structure prototypes it has not seen at least several hundred examples of in its training data. Hence Wrenformer may be unsuitable for discovering new prototypes.

<ProtoCountsWrenformerFailures />

## Element Prevalence vs Model Error

{#if mounted}
<ElementPrevalenceVsError />
{/if}

> @label:fig:element-prevalence-vs-error The y-axis is the average error of each model on each element averaged over all structures in the test set containing said element weighted by that structure's composition. The x-axis is each element's prevalence in the MP training set. We don't observe much correlation between element prevalence and model error. Instead, oxygen, the element with highest representation, is among the higher error elements for most models. Similarly, fluorine, the element with highest average error across models has very good training support at ~12 k samples, suggesting chemistry more than training coverage determines element errors.

## Hull Distance Box plot

{#if mounted}
<BoxHullDistErrors />
{/if}

> @label:fig:box-hull-dist-errors Box plot of interquartile ranges (IQR) of hull distance errors for each model. The whiskers extend to the 5th and 95th percentiles. The horizontal line inside the box shows the median. BOWSR has the largest median error while Voronoi RF has the largest IQR. Note that MEGNet and CGCNN are the only models with positive median. Their hull distance errors are biased towards more frequently predicting thermodynamic instability, explaining why they are closest to getting the overall number of stable structures in the test set right ([see rolling classification precision/recall plots](/preprint#fig:cumulative-precision-recall))
