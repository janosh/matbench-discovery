<script lang="ts">
  import MetricsTableMegnetCombos from '$figs/metrics-table-megnet-combos.svelte'
  import MetricsTableFirst10k from '$figs/metrics-table-first-10k.svelte'
  import RunTimeBars from '$figs/model-run-times-bar.svelte'
  import RocModels from '$figs/roc-models.svelte'
  import { browser } from '$app/environment'
  import MPRefEnergies from '$figs/mp-elemental-ref-energies.svelte'
  import WrenformerRollingMaeBatches from '$figs/rolling-mae-vs-hull-dist-wbm-batches-wrenformer.svelte'
  import CHGNetRollingMaeBatches from '$figs/rolling-mae-vs-hull-dist-wbm-batches-chgnet.svelte'
  import M3gnetRollingMaeBatches from '$figs/rolling-mae-vs-hull-dist-wbm-batches-m3gnet.svelte'
  import MegnetRollingMaeBatches from '$figs/rolling-mae-vs-hull-dist-wbm-batches-megnet.svelte'
  import CgcnnRollingMaeBatches from '$figs/rolling-mae-vs-hull-dist-wbm-batches-cgcnn.svelte'
  import VoronoiRfRollingMaeBatches from '$figs/rolling-mae-vs-hull-dist-wbm-batches-voronoi-rf.svelte'
  import HistClfPredHullDistModels from '$figs/hist-clf-pred-hull-dist-models.svelte'
  import SpacegroupSunburstWbm from '$figs/spacegroup-sunburst-wbm.svelte'
  import SpacegroupSunburstWrenformerFailures from '$figs/spacegroup-sunburst-wrenformer-failures.svelte'
  import ScatterLargestErrorsModelsMeanVsEachTrue from '$figs/scatter-largest-errors-models-mean-vs-each-true.svelte'
  import EAboveHullScatterWrenformerFailures from '$figs/e-above-hull-scatter-wrenformer-failures.svelte'
  import ProtoCountsWrenformerFailures from '$figs/proto-counts-wrenformer-failures.svelte'
  import { onMount } from 'svelte'

  let mounted = false
  onMount(() => (mounted = true))
</script>

# Supplementary Information

## Metrics for 10k materials predicted most stable

<MetricsTableFirst10k />

> @label:fig:metrics-table-first-10k An actual discovery campaign is unlikely to validate all of a model's stable predictions like we did in the [metrics table](/preprint#fig:metrics-table). Presumably it will rank model predictions from most to least stable and go down that list as far time and compute budgets permit. Assuming that increases in compute resources will allow average future discovery campaigns to grow in scope, we believe 10 k model validations to be a reasonable cutoff. To simulate this scenario, we calculated classification and regression metrics for the 10 k test set materials predicted to be most stable by each model.<br>
> We again show dummy performance in the bottom row. Note that each model is now evaluated on a different slice of the data, but this is still dummy performance across the whole dataset. CHGNet and M3GNet achieve a very impressive 83% and 80% precision, respectively. In concrete terms, this means in a discovery campaign that validates 10 k model predictions from a search pool of 257 k crystals which are chemically dissimilar from the training set and of which 16.7 % are stable, CHGNet and M3GNet would deliver 4 stable structures for every 5 predictions validated. In light of the significant resulting increase in stability hit rate, these models are well worth integrating into future materials searches.

## ROC Curves

{#if mounted}
<RocModels />
{/if}

> @label:fig:roc-models Receiver operating characteristic (ROC) curve for each model. TPR/FPR = true/false positive rate. FPR on the $x$-axis is the fraction of unstable structures classified as stable. TPR on the $y$-axis is the fraction of stable structures classified as stable. Points are colored by stability threshold $t$ which sweeps from $-0.4 \ \frac{\text{eV}}{\text{atom}} \leq t \leq 0.4 \ \frac{\text{eV}}{\text{atom}}$ above the hull. A material is classified as stable if the predicted E<sub>above hull</sub> lies below the stability threshold. Since all models predict E<sub>form</sub> (and M3GNet predicted energies are converted to formation energy before stability classification), they are insensitive to changes in the threshold $t$. M3GNet wins in area under curve (AUC) with 0.87, coming in 34% higher than the worst model Voronoi Random Forest. The diagonal 'No skill' line shows performance of a dummy model that randomly ranks material stability.

## Model Run Times

{#if mounted}
<RunTimeBars style="margin: 1em;" />
{/if}

> @label:fig:model-run-times-pie Creating this benchmark (excluding debugging runs) used a total of 3210 hours of compute time (mix of CPU and GPU, mostly CPU). Notably, the vast majority of that (2705 h) was used in the Bayesian optimization step of BOWSR + MEGnet.<br>
> Some bars have two sections. The bottom shows training time, the upper test time. If there's only one section (everything right of BOWSR + MEGnet) the model was not re-trained for this benchmark and the bar shows only the test time.

## Formation Energy MAE = Hull Distance MAE

To avoid potential confusion for people reading the code, we may in places calculate the formation energy MAE and report it as the MAE for the energy above the convex hull prediction. The former is more easily calculated but the two quantities are the same. The formation energy of a material is the difference in energy between a material and its constituent elements in their standard states. The distance to the convex hull is defined as the difference between a material's formation energy and the minimum formation energy of all possible stable materials made from the same elements. Since the formation energy of a material is used to calculate the distance to the convex hull, the error of a formation energy prediction directly determines the error in the distance to the convex hull prediction.

A further point of clarification: whenever we say distance to the convex hull we mean the signed distance that is positive for thermodynamically unstable materials above and negative for stable materials below the hull.

## MP Elemental Reference Energies

{#if mounted}
<MPRefEnergies />
{/if}

> @label:fig:mp-elemental-reference-energies WBM formation energies were calculated w.r.t. these Materials Project elemental reference energies ([queried on 2022-09-19](https://github.com/janosh/matbench-discovery/blob/main/data/mp/2022-09-19-mp-elemental-reference-entries.json)). Marker size indicates the number of atoms in the reference structure. Hover points for details.

## Classification Histograms using Model-Predicted Energies

{#if mounted}
<HistClfPredHullDistModels />
{/if}

> @label:fig:hist-clf-pred-hull-dist-models Similar to the [true hull distance histograms](/preprint#fig:hist-clf-true-hull-dist-models), these histograms show model stability classification as a function of the distance to the convex hull. The difference here being the $x$ axis shows model-predicted rather than DFT hull distances. Intuitively, it shows how often models misclassify as a function of how far they think a material is from the convex hull.

Note the CGCNN+P histogram is more strongly peaked than CGCNN's which agrees better with the actual DFT ground truth [distribution of hull distances](/about-the-data#--target-distribution) in our test set. This explains why CGCNN+P performs better as a regressor but also reveals how it can simultaneously perform worse as a classifier. By moving predictions closer to the stability threshold at 0 eV/atom above the hull, even small errors are significant enough to tip a classification over the threshold.

## WBM Batch Robustness as a Measure of Extrapolation Prowess

The figures below show the rolling MAE as a function of distance to the convex hull for each of the 5 WBM batches of elemental substitution. As a reminder, the WBM test set was generated in 5 successive batches, each applying another element replacement according to ISCD-mined chemical similarity scores to MP source structures and new stable crystals generated in the previous round. Naively, you'd expect model performance to degrade with increasing batch count as repeated substitutions should on average 'diffuse' deeper into uncharted regions of material space, requiring the model to extrapolate more. We observe this effect for some models but not others.

{#if mounted}

<div style="display: grid; grid-template-columns: 1fr 1fr; margin: 0 -1em 0 -4em;">
  <M3gnetRollingMaeBatches style="margin: -2em 0 0; height: 400px;" />
  <CHGNetRollingMaeBatches style="margin: -2em 0 0; height: 400px;" />
  <WrenformerRollingMaeBatches style="margin: -2em 0 0; height: 400px;" />
  <MegnetRollingMaeBatches style="margin: -2em 0 0; height: 400px;" />
  <VoronoiRfRollingMaeBatches style="margin: -2em 0 0; height: 400px;" />
  <CgcnnRollingMaeBatches style="margin: -2em 0 0; height: 400px;" />
</div>
{/if}

> @label:fig:rolling-mae-vs-hull-dist-wbm-batches-models On WBM batch 1, Wrenformer performs comparably to the SOTA UIPs M3GNet and CHGNet. However, the M3GNet and CHGNet UIPs show stronger extrapolative performance as they barely deteriorate in performance on later batches that move further away from the original training distribution. Wrenformer by contrast exhibits a pronounced increase in MAE with batch count.

MEGNet and CGCNN both incur a higher rolling MAE than Wrenformer across all 5 batches and across most or all of the hull distance range visible in these plots. However, similar to the UIPs, MEGNet and CGCNN exhibit very little degradation for higher batch counts. The fact that higher errors for later batches is specific to Wrenformer and Voronoi RF suggests that training only on composition and coarse-grained structural features (spacegroup and Wyckoff positions in the case of Wrenformer; coordination numbers, local environment properties, etc. in the case of Voronoi RF) alone is insufficient to learn an extrapolatable map of the PES.

Given its strong performance on batch 1, it is possible that given sufficiently diverse training data, Wrenformer could become similarly accurate to the UIPs across the whole PES landscape at substantially less training and inference cost. However, the loss of predictive forces and stresses may make Wrenformer unattractive for certain applications even then.

## Largest Errors vs DFT Hull Distance

{#if mounted}
<ScatterLargestErrorsModelsMeanVsEachTrue />
{/if}

> @label:fig:scatter-largest-errors-models-mean-vs-each-true The 200 structures with largest error averaged over all models vs their DFT hull distance colored by model disagreement (as measured by standard deviation in hull distance predictions from different models) and sized by number of training structures containing the least prevalent element (e.g. if a scatter point had composition FeO, MP has 6.6k structures containing Fe and 82k containing O so its size would be set to 6.6k). Thus smaller points have less training support. This plot suggests all models are biased to predict low energy and perhaps fail to capture certain physics resulting in highly unstable structures. This is unsurprising considering MP training data mainly consists of low energy structures.<br>
> It is also possible that some of the blue points with large error yet good agreement among models are in fact accurate ML predictions for a DFT relaxation gone wrong.

## MEGNet formation energies from UIP-relaxed structures

{#if mounted}
<MetricsTableMegnetCombos select={[`model`, `MEGNet`, `CHGNet`, `M3GNet`, `CHGNet + MEGNet`, `M3GNet + MEGNet`]} />
{/if}

> @label:fig:metrics-table-megnet-combos This table shows metrics obtained by combining MEGNet with both UIPs. The metrics in rows labeled M3GNet + MEGNet and CHGNet + MEGNet are the result of passing M3GNet/CHGNet-relaxed structures into MEGNet for formation energy prediction. Both combos perform worse than using the respective UIPs on their own with a more pronounced performance drop from CHGNet to CHGNet + MEGNet than M3GNet to M3GNet + MEGnet. This suggests MEGNet has learned no additional knowledge of the PES that is not already present in the UIPs. However, both combos perform better than MEGNet on its own, demonstrating that UIP relaxation provides real utility at very low cost for any downstream structure-dependent analysis.

The UIPs M3GNet and CHGNet are both trained to predict DFT energies (including/excluding MP2020 energy corrections for CHGNet/M3GNet) while MEGNet is trained to predict formation energies.

The two types of target energies should in principle be equivalent since the difference between raw DFT energies and (corrected) formation energies is just a linear transformation of subtracting the elemental reference energies weighted by composition. In that sense, changing the target from DFT energies to formation energies merely constitutes a change of gauge. But in practice, one might expect the problem of stability prediction to be easier when trained on formation energies as the model output is one less step removed from the ultimate target of crystal stability. Empirical evidence so far suggests the opposite, however. The UIP + MEGNet combos both perform worse than UIPs by themselves, suggesting using formation energy as targets offers no benefits.

We highlight this here to refute the suggestion that training on raw DFT energies results in poorly calibrated predictions for deriving formation energies and subsequent stability predictions.

## Spacegroup prevalence in Wrenformer failure cases

{#if mounted}

<div style="display: flex; gap: 1em; justify-content: space-around;">
<SpacegroupSunburstWrenformerFailures />
<SpacegroupSunburstWbm />
</div>
<EAboveHullScatterWrenformerFailures style="height: 300; width: 300;" />
{/if}

> @label:fig:spacegroup-prevalence-wrenformer-failures The left spacegroup sunburst shows spacegroup 71 is by far the dominant lattice symmetry among the 941 Wrenformer failure cases where $E_\text{above hull,DFT} < 1$ and $E_\text{above hull,Wrenformer} > 1$ (points inside the shaded rectangle). On the right side for comparison is the spacegroup sunburst for the entire WBM test set.

Looking at the occurrence counts of isopointal prototypes in the shaded rectangle and comparing them with the occurrence of those same prototypes in the MP training data counts, we find almost no support for failing structure prototypes. This suggests the reason Wrenformer fails so spectacularly on these structures is that it cannot deal with structure prototypes it has not seen at least several hundred examples of in its training data. Hence Wrenformer may be unsuitable for discovering new prototypes.

<ProtoCountsWrenformerFailures />
