<script>
  import { browser } from '$app/environment'
  import { onMount } from 'svelte'
  import { repository as repo } from '$site/package.json'
  import BoxHullDistErrors from '$figs/box-hull-dist-errors.svelte'
  import CgcnnRollingMaeBatches from '$figs/rolling-mae-vs-hull-dist-wbm-batches-cgcnn.svelte'
  import CHGNetRollingMaeBatches from '$figs/rolling-mae-vs-hull-dist-wbm-batches-chgnet.svelte'
  import CumulativeMae from '$figs/cumulative-mae.svelte'
  import CumulativePrecisionRecall from '$figs/cumulative-precision-recall.svelte'
  import EachScatterModels from '$figs/each-scatter-models-5x2.svelte'
  import ElementPrevalenceVsError from '$figs/element-prevalence-vs-error.svelte'
  import HistClfPredHullDistModels from '$figs/hist-clf-pred-hull-dist-models-5x2.svelte'
  import HullDistScatterWrenformerFailures from '$figs/hull-dist-scatter-wrenformer-failures.svelte'
  import LargestErrorScatterSelect from './largest-error-scatter-select.svelte'
  import M3gnetRollingMaeBatches from '$figs/rolling-mae-vs-hull-dist-wbm-batches-m3gnet.svelte'
  import MegnetRollingMaeBatches from '$figs/rolling-mae-vs-hull-dist-wbm-batches-megnet.svelte'
  import MetricsTable from '$figs/metrics-table.svelte'
  import MetricsTableFirst10k from '$figs/metrics-table-first-10k.svelte'
  import MetricsTableMegnetUipCombos from '$figs/metrics-table-uip-megnet-combos.svelte'
  import ProtoCountsWrenformerFailures from '$figs/proto-counts-wrenformer-failures.svelte'
  import RocModels from '$figs/roc-models-all-in-one.svelte'
  import RollingMaeVsHullDistModels from '$figs/rolling-mae-vs-hull-dist-models.svelte'
  import RunTimeBars from '$figs/model-run-times-bar.svelte'
  import SpacegroupSunburstWbm from '$figs/spacegroup-sunburst-wbm.svelte'
  import SpacegroupSunburstWrenformerFailures from '$figs/spacegroup-sunburst-wrenformer-failures.svelte'
  import VoronoiRfRollingMaeBatches from '$figs/rolling-mae-vs-hull-dist-wbm-batches-voronoi-rf.svelte'
  import WrenformerRollingMaeBatches from '$figs/rolling-mae-vs-hull-dist-wbm-batches-wrenformer.svelte'

  let mounted = false
  onMount(() => (mounted = true))
</script>

## Introduction

Material science can be viewed as a combinatorial problem of mixing and arranging different atoms to leverage the complex range of properties that emerge.
To date, $\sim 10^5$ combinations have been tested experimentally
[@bergerhoff_inorganic_1983; @belsky_new_2002] and $\sim 10^7$ have been simulated
[@jain_commentary_2013; @saal_materials_2013; @curtarolo_aflow_2012; @draxl_nomad_2018].
@davies_computational_2016 identified $\sim 10^{10}$ possible quaternary materials allowed by electronegativity and charge-balancing rules.
The space of quinternaries and higher is even less explored, leaving vast numbers of potentially useful materials to be discovered.
The discovery of new materials is a key driver of technological progress and lies on the path to more efficient solar cells, lighter and longer-lived batteries, smaller and more efficient transistor gates just to name a few.
In light of global warming, these advances cannot come fast enough. Any speed-up new methods might yield should be leveraged to their fullest extent.

Despite significant advances in empirical, theoretical and computational materials science, discovering new materials still requires complex calculations, labor-intensive trial-and-error experimentation, and often happens fortuitously rather than through rational design.
Machine learning (ML) methods efficiently extract and distill trends from huge datasets, can handle high dimensionality, multiple objectives, uncertainty, and noisy or sparse data, making them powerful additions to the computational materials science tool set.

ML models are less accurate and reliable but orders of magnitude faster than ab-initio simulation.
This makes them most suitable for use in high-throughput (HT) searches to triage more expensive, higher-fidelity simulation methods.
The use of neural networks for learning the Kohn-Sham density-functional theory (DFT) potential energy surface (PES) can be traced as far back as [@behler_generalized_2007].
This work kicked off rapid advances and significant efforts to fit ever more sophisticated ML models to known samples of the PES.
Initially, most of these models were trained and deployed as interatomic potentials to study known materials of interest, a workflow that requires curating custom training data for each new system of interest [@bartok_machine_2018; @deringer_general-purpose_2020].
As larger and more diverse datasets have emerged from initiatives like the Materials Project (MP) [@jain_commentary_2013], AFLOW [@curtarolo_aflow_2012] or the Open Quantum Materials Database (OQMD) [@saal_materials_2013], researchers have begun to train so-called universal models that cover 90 or more of the most-application relevant elements in the periodic table.
This opens up the prospect of ML-guided materials discovery to increase the hit rate of stable crystals and speed up DFT- and expert-driven searches.

Progress in ML for materials is often measured according to performance on standard benchmark data sets.
As ML models have grown in complexity and applicability, benchmark datasets need to grow with them to accurately measure their usefulness.
However, due the the rapid pace of the field and the variety of possible approaches for framing the discovery problem, no large-scale benchmark yet exists for measuring the ability of ML to accelerate materials discovery.
As a result, it is unclear which methodologies or models are best suited for this task.
Recent works focusing on prospective computational materials discovery have proposed strategies based on one-shot coordinate free predictors [@goodall_rapid_2022], iterative Bayesian optimizers [@zuo_accelerating_2021], and universal interatomic potentials (UIP) [@chen_universal_2022; @deng_chgnet_2023; @batatia_mace_2023].
These papers deploy their respective model on specific systems and custom datasets to highlight certain strengths but have yet to be compared in a systematic standardized manner that allows them to be ranked by their ability to accelerate materials discovery.
Our work aims to identify the state-of-the-art (SOTA) model by proposing a novel evaluation framework that closely simulates a real-world discovery campaign guided by ML models.

Specifically, we designed a benchmark task that tackles three central challenges. We believe these challenges to be requirements when seeking to justify experimental validation of ML predictions:

1. **Mimic Real Applications**: Idealized and overly simplified benchmarks can fail to reflect the real-world challenges a model faces when used in an actual discovery campaign.
   This can result in a disconnect between benchmark metrics and real-world performance.
   Possible reasons for this include choosing the wrong target [@bartel_critical_2020] or picking an unrepresentative train/test split [@wu_moleculenet_2018; @kpanou_robustness_2021]. In the case of materials discovery, formation energies - although widely used as regression targets - do not indicate thermodynamic stability.
   That is determined by the distance to the convex hull spanned by competing phases in the same chemical system.
   Moreover, ML models relying on relaxed crystal structures as input render any discovery pipeline circular since obtaining relaxed structures requires computationally expensive DFT simulations, thereby depending on the very process we intend to accelerate.

2. **Opportunity cost**: Accurate regressors are susceptible to unexpectedly high false-positive rates if those accurate predictions lie close to the decision boundary.
   Looking purely at global metrics like $\text{MAE}$, $\text{RMSE}$ and $R^2$ can give practitioners a false sense of security about their model's reliability.
   Failed experiments incur a high opportunity cost by wasting lab time and resources.

3. **Scalability**: Future discovery efforts are likely to encounter large data regimes.
   Small benchmarks can lack chemical diversity, obfuscate poor scaling relations or poor out-of-distribution (OOD) performance.
   For instance, random forests achieve excellent performance on small datasets but are typically outperformed by neural networks on large datasets due to the benefits of representation learning [@goodall_predicting_2020].

Two prior efforts worth highlighting that have partially addressed the above challenges are Matbench [@dunn_benchmarking_2020] and the Open Catalyst Project (OCP) [@chanussot_open_2021].

By providing a standardized collection of 13 datasets ranging in size from ~300 to ~132,000 samples from both DFT and experimental sources, Matbench addresses the scalability challenge, highlighting how model performance changes as a function of data regime.
Matbench helped focus the field of ML for materials, increase comparability across papers and provide a quantitative measure of progress in the field.
Importantly, all tasks were exclusively concerned with the properties of known materials.
We believe a task that simulates a materials discovery campaign by requiring materials stability prediction from unrelaxed structures to be a missing piece here.

OCP is a large-scale initiative to discover substrate-adsorbate combinations that catalyze key industrial reactions that process said adsorbates into more useful products.
The OCP has released two data sets thus far, OCP20 [@chanussot_open_2021] and OCP22 [@tran_open_2022], for training and benchmarking ML models.
OCP certainly addressed challenge 1 of closely mimicking a real-world problem.
They have recently shown that despite not reaching their target accuracy to entirely replace DFT, using ML in conjunction with confirmatory DFT calculations dramatically speeds up their combinatorial screening workflow [@lan_adsorbml_2023].

We believe that addressing these three challenges will enable future ML-guided discovery efforts to expand materials databases to confidently select appropriate models and methodology.
Our initial findings show universal potentials to outperform all other methodologies we tested both in accuracy and robustness.

### Evaluation Framework for Materials Discovery

We propose a novel evaluation framework that places no constraints on the type of training data models use.
We only enforce that at test time, all models must make predictions on the convex hull distance of the relaxed structure with only the unrelaxed structure as input.
This setup avoids circularity in the discovery process, as unrelaxed structures can be cheaply enumerated through elemental substitution methodologies and do not contain information inaccessible in a prospective discovery campaign.
We choose to focus on the relaxed structure's convex hull distance as a measure of thermodynamic stability rather than the formation energy as it informs the decision on whether to pursue a potential candidate crystal.
This decision was also motivated by @bartel_critical_2020 who found even composition-only models capable of predicting DFT formation energies with useful accuracy.
However, when tasking those same models with predicting decomposition enthalpy, performance deteriorated sharply.
This insight meant that ML models are much less useful than DFT for discovering new inorganic solids.
Moreover, the qualitative leap in performance from Roost [@goodall_predicting_2020], the best compositional model benchmarked, to CGCNN [@xie_crystal_2018], the single structural model they tested, highlights that structure plays a crucial role in determining material stability.
Hence our decision to restrict test-time model input to unrelaxed structures to avoid the above-mentioned circularity and measure true prospective utility while still allowing for subtle atomic-configuration-informed energy estimates.
We recognize that our ground truth is ignorant of entropic stabilization and metastable states.
While these factors influence experimental synthesizability, we have to neglect them when simulating an HT search since the convex hull distance as predicted by zero-temperature DFT is the best proxy of true crystal stability available for large datasets.

For the first realization of this framework we use the Materials Project (MP) [@jain_commentary_2013] database release ([v2022.10.28](https://docs.materialsproject.org/changes/database-versions) as our training set, and the unrelaxed structures in the WBM dataset [@wang_predicting_2021] as our test set.

#### Materials Project Training Set

The Materials Project is a well-known effort to calculate the properties of all inorganic materials using high-throughput ab-initio methods.
Seeded from a subset of the Inorganic Crystal Structure Database (ICSD) [@allen_crystallographic_1999], the initial release of the database consisted of ~9 k crystals.
At the time of writing, the Materials Project database [@jain_commentary_2013] has grown to ~154 k crystals, covering diverse chemistries and providing relaxed and initial structures as well as the relaxation trajectory for every entry.

Our benchmark defines the training set as all data available from the [2022.10.28 MP release](https://docs.materialsproject.org/changes/database-versions#v2022.10.28).
We recorded a snapshot of energies, forces, stresses and magnetic moments for all MP ionic steps on 2023-03-15 as the canonical training set for v1 of Matbench Discovery, and provide convenience functions through our [Python package](https://pypi.org/project/matbench-discovery) for easily feeding that data into future model submissions to our benchmark.
This flexibility is intended to allow authors to experiment with and exploit the large variety of available data.

#### WBM Test Set

The WBM data set [@wang_predicting_2021] consists of ~257 k structures generated via chemical similarity-based elemental substitution of MP source structures followed by DFT relaxation and calculating each crystal's convex hull distance.
Which element replaces an existing one in a given source structure was determined by random sampling according to the weights in a chemical similarity matrix data-mined from the ICSD [@glawe_optimal_2016].

The WBM authors performed 5 iterations of this substitution process (we refer to these as batches).
After each step, the newly generated structures found to be thermodynamically stable after DFT relaxation flow back into the source pool to partake in the next round of substitution.
This split of the data into batches of increasing substitution count is a unique and compelling feature of the test set as it allows out-of-distribution (OOD) testing by seeing if model performance degrades with substitution count.
A higher number of elemental substitutions on average carries the structure further away from the region of material space covered by the MP training set (see @fig:rolling-mae-vs-hull-dist-wbm-batches-models for details).
Whilst this batch information makes the WBM data set an exceptionally useful data source for examining the extrapolation performance of ML models, our evaluation primarily looks at metrics that consider all batches as a single test set.

Throughout this work, we define stability as being on or below the convex hull of the MP training set, ~42k out of ~257k materials in WBM satisfy this criterion.
Our code treats the stability threshold as a dynamic parameter for future more detailed performance analysis at different thresholds.
For initial analysis in this direction, see @fig:roc-models in the SI.

As WBM explores regions of materials space not well sampled by MP, many of the discovered materials that lie below MP's convex hull are not stable relative to each other.
Of the ~42k that lie below the MP convex hull less than half, or around ~20k, remain on the convex hull when merging the MP and WBM hulls.
This observation suggests that many WBM structures are repeated samples in the same new chemical spaces.
It also highlights a critical aspect of this benchmark in that we purposely operate on an incomplete convex hull.
Only current knowledge is accessible to a real discovery campaign and our metrics are designed to reflect this.

### Models

To test a wide variety of methodologies proposed for learning the potential energy landscape, our initial benchmark release includes 10 models.

1. **CHGNet** [@deng_chgnet_2023] (UIP-GNN) - CHGNet is a UIP for charge-informed atomistic modeling.
   Its distinguishing feature is that it was trained to predict magnetic moments on top of energies, forces and stresses in the MPTrj dataset consisting of relaxation trajectories for ~1.5 million MP structures.
   By modeling magnetic moments, CHGNet learns to accurately represent the orbital occupancy of electrons which allows it to predict both atomic and electronic degrees of freedom.

2. **M3GNet** [@chen_universal_2022] (UIP-GNN) - M3GNet is a GNN-based universal (as in full periodic table) interatomic potential (UIP) for materials trained on up to 3-body interactions in the initial, middle and final frame of MP DFT relaxations.
   The model takes the unrelaxed input and emulates structure relaxation before predicting energy for the pseudo-relaxed structure.

3. **MACE** [@batatia_mace_2023] (UIP-GNN) - MACE builds upon the recent advances [@batzner_equivariant_2022; @thomas_tensor_2018] in equivariant neural network architectures by proposing an approach to computing high N-body order features in an efficient manner via Atomic Cluster Expansion [@drautz_atomic_2019].
   Unlike the other UIP models considered MACE was primarily developed for molecular dynamics of single material systems and not the universal use case studied here.
   It is the only equivariant model we tested.

4. **ALIGNN** [@choudhary_atomistic_2021] (GNN) - The Atomistic Line Graph Neural Network (ALIGNN) is a message passing GNN architecture that takes as input both the interatomic bond graph and a line graph corresponding to 3-body bond angles.
   The ALIGNN architecture involves a global pooling operation which means that it is ill-suited to force-field applications.
   To address this the ALIGNN-FF model was later introduced without global pooling [@choudhary_unified_2023].

5. **MEGNet** [@chen_graph_2019] (GNN) - MatErials Graph Network is another GNN similar to CGCNN for material properties of relaxed structures that also updates the edge and global features (like pressure, temperature, entropy) in its message passing operation.
   This work showed that learned element embeddings encode periodic chemical trends and can be transfer-learned from large data sets (formation energies) to predictions on small data properties (band gaps, elastic moduli).

6. **CGCNN** [@xie_crystal_2018] (GNN) - The Crystal Graph Convolutional Neural Network (CGCNN) was the first neural network model to directly learn 8 different DFT-computed material properties from a graph representing the atoms and bonds in a periodic crystal.
   CGCNN was among the first to show that just like in other areas of ML, given large enough training sets, neural networks can learn embeddings that outperform human-engineered structure features directly from the data.

7. **CGCNN+P** [@gibson_data-augmentation_2022] (GNN) - This work proposes simple, physically motivated structure perturbations to augment stock CGCNN's training data of relaxed structures with structures resembling unrelaxed ones but mapped to the same DFT final energy.
   Here we chose $P=5$, meaning the training set is augmented with 5 random perturbations of each relaxed MP structure mapped to the same target energy.
   In contrast to all other structure-based GNNs considered in this benchmark, CGCNN+P is not attempting to learn the Born-Oppenheimer potential energy surface.
   The model is instead taught the PES as a step-function that maps each valley to its local minimum.
   The idea is that during testing on unrelaxed structures, the model will predict the energy of the nearest basin in the PES.
   The authors confirm this by demonstrating a lowering of the energy error on unrelaxed structures.

8. **Wrenformer** (Transformer) - For this benchmark, we introduce Wrenformer which is a variation on the coordinate-free Wren model [@goodall_rapid_2022] constructed using standard QKV-self-attention blocks [@vaswani_attention_2017] in place of message-passing layers.
   This architectural adaptation reduces the memory usage allowing the architecture to scale to structures with greater than 16 Wyckoff positions.
   Like its predecessor, Wrenformer is a fast coordinate-free model aimed at accelerating screening campaigns where even the unrelaxed structure is a priori unknown.
   The key idea is that by training on the coordinate anonymized Wyckoff positions (symmetry-related positions in the crystal structure), the model learns to distinguish polymorphs while maintaining discrete and computationally enumerable inputs.
   The central methodological benefit of an enumerable input is that it allows users to predict the energy of all possible combinations of spacegroup and Wyckoff positions for a given composition and maximum unit cell size.
   The lowest-ranked prototypes can then be fed into downstream analysis or modeling.

9. **BOWSR** [@zuo_accelerating_2021] (BO-GNN) - BOWSR combines a symmetry-constrained Bayesian optimizer (BO) with a surrogate energy model to perform an iterative exploration-exploitation-based search of the potential energy landscape.
   Here we use MEGNet [@chen_graph_2019] for the energy model as proposed in the original work. The high sample count needed to explore the PES with BO makes this by far the most expensive model tested.

10. **Voronoi RF** [@ward_including_2017] (Fingerprint) - A random forest trained to map a combination of composition-based Magpie features [@ward_general-purpose_2016] and structure-based relaxation-robust Voronoi tessellation features (effective coordination numbers, structural heterogeneity, local environment properties, ...) to DFT formation energies.
    This fingerprint-based model predates most deep learning for materials but significantly improved over earlier fingerprint-based methods such as the Coulomb matrix [@rupp_fast_2012] and partial radial distribution function features [@schutt_how_2014].
    It serves as a baseline model to see how much value the learned featurization of deep learning models can extract from the increasingly large corpus of available training data.

## Results

{#if mounted}
<MetricsTable />
{/if}

> @label:fig:metrics-table Classification and regression metrics for all models tested on our benchmark ranked by F1 score.
> The heat map ranges from yellow (best) to blue (worst) performance.
> DAF = discovery acceleration factor (see text), TPR = true positive rate, TNR = false negative rate, MAE = mean absolute error, RMSE = root mean squared error.
> The dummy classifier uses the 'scikit-learn' 'stratified' strategy of randomly assigning stable/unstable labels according to the training set prevalence.
> The dummy regression metrics MAE, RMSE and $R^2$ are attained by always predicting the test set mean.
> The Voronoi RF, CGCNN and MEGNet models are seen to be worse than the dummy result on regression metrics but better on some of the classification metrics, highlighting the importance of looking at the right metrics for the task at hand to gauge model performance.
>
> <details>
> <summary>Table glossary</summary>
>
> - DAF = discovery acceleration factor
> - TPR = true positive rate, the fraction of stable structures correctly predicted as stable
> - TNR = true negative rate, the fraction of unstable structures correctly predicted as unstable
> - MAE = mean absolute error
> - RMSE = root mean squared error
> - GNN = graph neural network
> - UIP = universal interatomic potential
> - BO = Bayesian optimization
> - RF = random forest
> - +P = training data augmentation using random structure perturbations
>
> </details>

@Fig:metrics-table shows performance metrics for all models included in the initial release of Matbench Discovery.
CHGNet takes the top spot on all metrics except true positive rate (TPR) and emerges as the current SOTA for ML-guided materials discovery.
The discovery acceleration factor (DAF) measures how many more stable structures a model found compared to the dummy discovery rate of 43k / 257k $\approx$ 16.7% achieved by randomly selecting test set crystals.
The maximum possible DAF is the inverse of the dummy discovery rate which on our dataset is ~6.
The current SOTA of 3.06 achieved by CHGNet leaves room for improvement.
However, evaluating models on their 10k most stable predictions only as shown in @fig:metrics-table-first-10k, CHGNet achieves an impressive DAF of 4.96 approaching optimal performance.

We find a large performance gap between models that make one-shot predictions directly from unrelaxed inputs (MEGNet, Wrenformer, CGCNN, CGCNN+P, ALIGNN, Voronoi RF) compared to UIPs that predict forces and use them to emulate DFT relaxation (CHGNet, M3GNet, MACE).
While the F1 scores and DAFs of force-free models are less unaffected, their coefficients of determination ($R^2$) are significantly worse.
Of the force-free models, only ALIGNN, BOWSR and CGCNN+P achieve positive $R^2$.
Negative $R^2$ means model predictions explain the observed variation in the data less than simply predicting the test set mean.
In other words, these models are not predictive in a global sense (across the full dataset range).
However, even models with negative $R^2$ can be locally good in the positive and negative tails of the test set hull distance distribution.
They suffer most in the mode of our hull distance distribution near the stability threshold of 0 eV/atom above the hull.
This reveals an important shortcoming of $R^2$ as a metric for classification tasks like stability prediction.

The results for M3GNet and MACE depart from the trend that F1 is rank-correlated with the other classification metrics.
Of all models, M3GNet achieves the highest true positive rate (TPR) but an unusually low true negative rate (TNR).
A similar trend is seen for MACE. @fig:rolling-mae-vs-hull-dist-models provides a visual understanding of this observation.
M3GNet and MACE have the lowest rolling mean of the absolute errors (rolling MAE) as a function of hull distance for materials above the convex hull (see right half of plot) but incur comparably large errors for materials below the hull (left half of plot).
Since $\text{TPR} = \frac{\text{TP}}{\text{TP} + \text{FN}}$, lower error for energies above the hull increases both TN and decreases FP, resulting in the high TPR values observed.

The reason CGCNN+P achieves better regression metrics than CGCNN but is still worse as a classifier becomes apparent from @fig:hist-clf-pred-hull-dist-models by noting that the CGCNN+P histogram is more sharply peaked at the 0 hull distance stability threshold.
This causes even small errors in the predicted convex hull distance to be large enough to invert a classification.
Again, this is evidence to choose carefully which metrics to optimize.
Regression metrics are far more prevalent when evaluating energy predictions. However, our benchmark treats energy predictions as merely means to an end to classify compound stability.
Improvements in regression accuracy are of limited use to materials discovery in their own right unless they also improve classification accuracy.
Our results demonstrate that this is not a given.

{#if mounted}
<CumulativePrecisionRecall />
{/if}

> @label:fig:cumulative-precision-recall CHGNet achieves the highest cumulative precision and recall at any point during a simulated discovery campaign.
> A discovery campaign consists of using a model to rank test set materials by their predicted energy above the training set convex hull.
> A higher fraction of correct stable predictions corresponds to higher precision and fewer stable materials overlooked correspond to higher recall.
> This figure highlights how different models perform better or worse depending on the length of the discovery campaign.
> The UIP models (CHGNet, M3GNet, MACE) are seen to offer significantly improved precision on shorter campaigns as they are less prone to early false positive predictions.

@Fig:cumulative-precision-recall has models rank materials by model-predicted hull distance from most to least stable; materials furthest below the known hull at the top, materials right on the hull at the bottom.
For each model, we iterate through that list and calculate at each step the precision and recall of correctly identified stable materials.
This simulates exactly how these models would be used in a prospective materials discovery campaign and reveals how a model's performance changes as a function of the discovery campaign length. As a practitioner, you have a certain amount of resources available to validate model predictions. These curves allow you to read off the best model given these conditions and based on the optimal trade-off between fewer false positives (precision) or fewer negatives (recall) for the discovery task at hand.
In this case, it so happens that CHGNet achieves the highest precision _and_ recall at any number of screened materials.

A line terminates when a model believes there are no more materials in the WBM test set below the MP convex hull.
The dashed vertical line shows the actual number of stable structures in our test set.
All models are biased towards stability to some degree as they all overestimate this number, most of all BOWSR by 133%, least of all MEGNet by 30%.
This is only a problem in practice for exhaustive discovery campaigns that validate _all_ stable predictions from a model.
More frequently, model predictions will be ranked most-to-least stable and validation stops after some pre-determined compute budget is spent, say, 10k DFT relaxations.
In that case, the concentration of false positive predictions that naturally accumulates near the less stable end of the candidate list can be avoided with no harm to the campaign's overall discovery rate (see @fig:metrics-table-first-10k a similar metrics table considering only the 10k materials predicted by each model to be furthest below the known convex hull).

The diagonal Optimal Recall line would be achieved if a model never made a false negative prediction and stopped predicting stable crystals exactly when the true number of stable materials is reached.
Zooming in on the top-left corner of the precision plot, we observe that CHGNet is the only model without a sudden drop in precision right at the start of the discovery campaign.
This means CHGNet is the only model whose first few hundred most stable materials are largely actually stable.
MACE in particular suffers from a large number of initial failure cases.
These are unstable materials whose energy MACE underpredicts by several eV/atom (see @fig:each-scatter-models), resulting in MACE's initial precision starting at 0 before quickly climbing to ~0.8.
We experienced this behavior with multiple independently trained MACE variants; even more so with a checkpoint we received directly from the MACE authors which was trained on the M3GNet training set. The results shown here are from a superior MACE we trained ourselves on the much larger MPtrj dataset.
Even M3GNet, the other universal potential and 2nd best model, exhibits the early-on precision drop.
As a result, CHGNet has a strong lead over M3GNet until reaching ~3k screened materials.
From there, CHGNet and M3GNet slowly converge until they almost tie at a precision of 0.52 after ~56k screened materials.
At that point, CHGNet's list of stable predictions is exhausted while M3GNet continues, dropping in precision to 0.45 at 76 k, attributable to many false positives near the end of the list of stable predictions.

All force-free models exhibit a much worse case of early-on precision drop, falling to 0.6 or less in the first 5k screened materials. Many of these models (all except BOWSR, Wrenformer and Voronoi RF) display an interesting hook shape in their cumulative precision, recovering again slightly in the middle of the simulated campaign between  5k and up to  30k before dropping again until the end.

{#if mounted}
<RollingMaeVsHullDistModels />
{/if}

> @label:fig:rolling-mae-vs-hull-dist-models Universal potentials are more reliable classifiers because they exit the red triangle earliest.
> These lines show the rolling MAE on the WBM test set as the energy to the convex hull of the MP training set is varied, lower is better.
> The large red 'triangle of peril' shows where the models are most likely to misclassify structures.
> As long as a model's rolling MAE remains inside the triangle, its mean error is larger than the distance to the convex hull.
> If the model's error for a given prediction happens to point towards the stability threshold at 0 eV from the hull (the plot's center), its average error will change the stability classification of a material from true positive/negative to false negative/positive.
> The width of the 'rolling window' box indicates the width over which errors hull distance prediction errors were averaged.

@Fig:rolling-mae-vs-hull-dist-models provides a visual representation of the reliability of different models based on the rolling mean absolute error (MAE) of model-predicted hull distances as a function of DFT distance to the Materials Project (MP) convex hull.
The red-shaded area, referred to as the 'triangle of peril', emphasizes the zone where the average model error surpasses the distance to the stability threshold at 0 eV.
As long as the rolling MAE remains within this triangle, the model is most susceptible to misclassifying structures.
Because the average error is larger than the distance to the classification threshold at 0, it is large enough to flip a correct classification into an incorrect one (if the error happens to point toward the stability threshold).

The sooner a model leaves the triangle on the left side, the less likely it is to incorrectly predict stable structures as unstable, thereby reducing false negatives.
On the right side, early exits result in a lower likelihood of predicting unstable structures as stable, leading to fewer false positives.

On the left side, CHGNet exits the triangle first. We can expect fewer false negatives from CHGNet than from other models.
On the right side, M3GNet is first to exit, at a remarkably low ~40 meV/atom error.
M3GNet is overtaken by Wrenformer for highly unstable structures towards the right end of the plot.
Therefore, both of these models are expected to produce fewer false positives.

Overall, this visualization underscores the importance of considering a model's average error relative to the distance to the stability threshold when evaluating its performance as a material stability predictor.
Low overall hull distance error can be a misleading metric if that error is not particularly low at the decision boundary.

The imbalance between the left and right half of the plot shows that models are more prone to false negative predictions, even for very stable materials far below the known hull, than predicting very unstable materials as stable.

As the convex hull becomes more thoroughly sampled by future discovery, the fraction of unknown stable structures decreases, naturally leading to less enriched future test sets.
This has several implications, including allowing for higher maximum DAFs but also skewing the hull distance distribution towards positive values (i.e. the right half of @fig:rolling-mae-vs-hull-dist-models) as there are fewer undersampled chemical spaces. Consequently, to have accurate stability classifications, model accuracy will need to improve concomitantly.
For @fig:rolling-mae-vs-hull-dist-models, this means models need to be much more accurate to exit the shaded triangle in the left half of the plot.

## Discussion

We have demonstrated the effectiveness of ML-based triage in HT materials discovery and posit that the benefits of including ML in discovery workflows now clearly outweigh the costs.
@fig:metrics-table shows in a realistic benchmark scenario that several models achieve a discovery acceleration greater than 2.5 across the whole dataset and up to 5 when considering only the 10k most stable predictions from each model (@fig:metrics-table-first-10k).
When starting this project, we were unsure which is the most promising ML methodology for HT discovery.
Our findings demonstrate a clear superiority in accuracy and extrapolation performance of UIPs like CHGNet, M3GNet and MACE.
Modeling forces enables these models to chart a path through atomic configuration space closer to the DFT-relaxed structure from where a more informed final energy prediction is possible.
The consistently linear log-log learning curves observed in related literature [@vonlilienfeld_retrospective_2020] suggest that further decreases in the error of UIPs can be readily unlocked with increased training data.

Despite its fast pace of progress, the field of universal ML potentials is still in its early stages.
To realize the full potential of these models, we will need to pivot significant resources to generate large quantities of higher-than-PBE fidelity training data.
That said, we believe our work (e.g. @fig:rolling-mae-vs-hull-dist-wbm-batches-models) and others [@lan_adsorbml_2023] have caught early glimpses that sufficiently large training sets allow us to narrow the performance gap between in-domain and out-of-domain performance.
This should further increase the discovery acceleration factor, resulting in a virtuous cycle between ML models accelerating the creation of new training data while becoming more accurate and even stronger accelerators after retraining.
Meanwhile, the setup cost and learning curve of ML potentials are steadily decreasing, likely with a much lower floor than DFT given the intricacies of HT ab-initio simulation.
Combined with their already strong performance, this suggests that UIPs combined with data accumulation over the next few years may unlock our way to a comprehensive yet cheap map of the PES that truly commoditizes HT materials discovery.

Despite this optimistic outlook, not only has the predictive power UIPs attained for crystal stability over the last 3 years since @bartel_critical_2020 sobering 2020 analysis yet to enter common knowledge, but many remaining shortcomings invite further research.
For example, all models are biased towards stability, overestimating the correct number of stable materials in our test set anywhere between 30% and 133%.
Moreover, even the best models in our analysis still make large numbers of false positive predictions, even for materials far above the convex hull.
This indicates the need to train future models on even more diverse datasets that are less biased towards ground state crystals than the MP training set we used.
No such dataset built with MP-compatible DFT settings exists to our knowledge, inviting future efforts in this direction.

Looking beyond mere thermodynamic stability prediction at zero Kelvin for the purpose of materials discovery, future materials science will also require understanding and predicting the properties under varying environmental conditions such as finite temperature and pressure.
This is where interatomic potentials can unlock further utility. Their force predictions enable insights into dynamical properties.
One open question of particular relevance is the extent to which such models may aid research into the computational prediction of synthesis pathways.
Many existing approaches for reaction pathway prediction involve the use of heuristic rules for dealing with the significant added complexity of metastability alongside traditional ground state ab-initio data [@mcdermott_graph-based_2021; @aykol_rational_2021; @wen_chemical_2023].
These algorithms stand to benefit massively from more efficient estimates of the reaction energy barriers that future UIPs may provide.
Provided they reach sufficient accuracy, this may enable timely calculation of approximate reaction barriers [@aykol_thermodynamic_2018] and open up a whole new field to high-throughput inquiry.

## Acknowledgments

J.R. acknowledges support from the German Academic Scholarship Foundation ([Studienstiftung](https://wikipedia.org/wiki/Studienstiftung)).
A.A.L. acknowledges support from the Royal Society.
A.J. and K.A.P. acknowledge the US Department of Energy, Office of Science, Office of Basic Energy Sciences, Materials Sciences and Engineering Division under contract no. DE-AC02-05-CH11231 (Materials Project program KC23MP).
This work used computational resources provided by the National Energy Research Scientific Computing Center (NERSC), a U.S. Department of Energy Office of Science User Facility operated under Contract No. DE-AC02-05CH11231.

Our profound gratitude extends to Hai-Chen Wang, Silvana Botti and Miguel A. L. Marques for their valuable contribution in crafting and freely sharing the WBM data set.

We would like to thank Jason Blake Gibson, Shyue Ping Ong, Chi Chen, Tian Xie, Bowen Deng, Peichen Zhong and Ekin Dogus Cubuk for helpful discussions.

We thank Rokas Elijosius for assisting in the initial implementation of Wrenformer.

## Author Contributions

Janosh Riebesell: Methodology, Software, Data Curation, Training and Testing Models, Formal Analysis. Rhys Goodall: Conceptualization, Software, Formal Analysis. Anubhav Jain: Supervision. Philipp Benner: Software, Training Models. Kristin Persson: Supervision. Alpha Lee: Supervision.

## Code availability

We welcome further model submissions to our GitHub repo <https://github.com/janosh/matbench-discovery>.

## Data availability

We chose the latest Materials Project (MP) [@jain_commentary_2013] database release ([v2022.10.28](https://docs.materialsproject.org/changes/database-versions) at time of writing) as the training set and the WBM dataset [@wang_predicting_2021] available at <https://figshare.com/articles/dataset/22715158> as the test set for this benchmark.
A snapshot of every ionic step including energies, forces, stresses and magnetic moments in the MP database is available at <https://figshare.com/articles/dataset/23713842>.

## Supplementary Information

### Metrics for 10k materials predicted most stable

An actual discovery campaign is unlikely to validate all stable predictions coming from a given model as we did in @fig:metrics-table.
Presumably, it will rank model predictions from most to least stable and follow that list as far as time and compute budget permits.
Assuming that increases in compute resources will allow average future discovery campaigns to grow in scope, we believe 10k model validations to be a reasonable scope for average campaigns.
This is what @fig:metrics-table-first-10k simulates by calculating classification and regression metrics for the 10k test set materials predicted to be most stable by each model.
We again show dummy performance in the bottom row. Note that each model is now evaluated on a different slice of the data. However, the bottom row still shows dummy performance across the whole dataset.
CHGNet and M3GNet achieve an impressive precision of 0.83 and 0.8, respectively.
In concrete terms, this means in a discovery campaign that validates 10 k model predictions from a search pool of 257 k crystals that are chemically dissimilar from the training set and of which 0.167 are stable, CHGNet and M3GNet would deliver 4 stable structures for every 5 predictions validated.
We emphasize that in light of the significant resulting increase in stability hit rate, these models are well worth integrating into future materials searches.

{#if mounted}
<MetricsTableFirst10k />
{/if}

> @label:fig:metrics-table-first-10k Model metrics for the 10k materials predicted to be most stable by each model.
> Each model predicts a different 10k subset of the test set as furthest below the known convex hull.
> CHGNet achieves an impressive F1 score of 0.91 and a discovery acceleration factor (DAF) of 4.96, both approaching the optimal values of 1 and 6, respectively.

### ROC Curves

A material is classified as stable if the predicted $E_\text{above hull}$ lies below the stability threshold. Since all models predict $E_\text{form}$ (and M3GNet predicted energies are converted to formation energy before stability classification), they are insensitive to changes in the threshold $t$. The receiver operating characteristic (ROC) curve for each model is plotted in @fig:roc-models. M3GNet wins in the area under curve (AUC) with 0.87, coming in 1.34x higher than the worst model Voronoi Random Forest. The diagonal 'No skill' line shows the performance of a dummy model that randomly ranks material stability.

{#if mounted}
<RocModels />
{/if}

> @label:fig:roc-models Receiver operating characteristic (ROC) curve for each model. TPR/FPR = true/false positive rate. FPR on the $x$-axis is the fraction of unstable structures classified as stable. TPR on the $y$-axis is the fraction of stable structures classified as stable. Points are colored by stability threshold $t$ which sweeps from $-0.4 \ \frac{\text{eV}}{\text{atom}} \leq t \leq 0.4 \ \frac{\text{eV}}{\text{atom}}$ above the hull.

### Parity Plots

{#if mounted}
<EachScatterModels />
{/if}

> @label:fig:each-scatter-models Parity plot for each model's energy above hull predictions (based on their formation energy predictions) vs DFT ground truth, color-coded by log density of points.

@Fig:each-scatter-models shows that all models do well for materials far below the convex hull (left side of the plot). Performance for materials far above the convex hull is more varied with occasional underpredictions of the energy of materials far above the convex hull (right side). All models suffer most in the mode of the distribution at $x = 0$.

Two models stand out as anomalous to the general trends.

Wrenformer is the only model with a large number of severe energy overpredictions at $x = 0$ going up the $y$-axis. We investigate these failure cases in more detail in @fig:wrenformer-failures and find these overpredictions to be dominated by spacegroup 71 with poor representation in the training data.

The other anomalous model is MACE with several severe underpredictions at $x = 0$ going down the $y$-axis. We investigated these points for common traits in composition or crystal symmetry but noticed no pattern.

Beyond these MACE outliers visible in the plot, MACE exhibits another rare but reproducible type of failure case in which the final predicted energy after relaxation is off by several orders of magnitude. The largest 'derailed' prediction was $-10^{22}$ eV/atom for `wbm-3-31970` (formula H$_2$Ir). In each case, the MACE relaxation exhausted the maximum number of ionic steps set to 500 and caused a volume implosion from initial cell volumes of hundreds to relaxed cell volumes of tens of cubic Angstrom. Using the checkpoint trained on the M3GNet dataset which we received from the MACE authors, this failure mode occurred for several hundred of the 250k test set crystals. Using the checkpoint we trained ourselves on the MPTrj dataset, it affects only 44 test crystals, suggesting that these holes in the MACE PES can perhaps be fully plugged by further increasing the training set or even changing the loss function.
Further analysis is ongoing.
Since these derailed values are easily identified in practice when actually performing a prospective discovery campaign, we excluded them from the MACE parity plat and all other downstream analyses.

### Model Run Times

{#if mounted}
<RunTimeBars style="margin: 1em;" />
{/if}

> @label:fig:model-run-times-pie Creating this benchmark (excluding debugging runs) used a total of 3210 hours of compute time (mix of CPU and GPU, mostly CPU). Notably, the vast majority of that (2705 h) was used in the Bayesian optimization step of BOWSR.<br>
> Some bars have two sections. The bottom shows training time, the upper test time. If there's only one section the model was not re-trained for this benchmark and the bar shows only the test time.

### Spacegroup prevalence in Wrenformer failure cases

{#if mounted}

<div style="display: flex; gap: 1em; justify-content: space-around; flex-wrap: wrap; margin-bottom: 2em;">
<SpacegroupSunburstWrenformerFailures />
<SpacegroupSunburstWbm />
</div>
<HullDistScatterWrenformerFailures />
{/if}

> @label:fig:wrenformer-failures Symmetry analysis of the 941 Wrenformer failure cases in the shaded rectangle defined by $E_\text{DFT hull dist} < 1$ and $E_\text{ML hull dist} > 1$. Sunburst plot of spacegroups shows that close to 80% of severe energy overestimations are orthorhombic with spacegroup 71. The table on the right shows occurrence counts of exact structure prototypes for each material in the sunburst plot as well as their corresponding prevalence in the training set.

@Fig:wrenformer-failures shows 456 + 194 ($\sim$ 70%) of the failure cases in the shaded rectangle are two prototypes in spacegroup 71.
The occurrence of those same prototypes in the MP training set shows almost no data support for the failing prototypes.
This suggests the reason Wrenformer fails so spectacularly on these structures is that it cannot deal with structure prototypes it has not seen at least several hundred examples of in its training data.
This suggests that there are stronger limitations on how much the discrete Wyckoff-based representation can extrapolate to new prototypes compared to the smooth local-environment-based inputs to GNN-type models.

<ProtoCountsWrenformerFailures />

### Hull Distance Box plot

{#if mounted}
<BoxHullDistErrors />
{/if}

> @label:fig:box-hull-dist-errors Box plot of interquartile ranges (IQR) of hull distance errors for each model. The whiskers extend to the 5th and 95th percentiles. The horizontal line inside the box shows the median. BOWSR has the largest median error, while Voronoi RF has the largest IQR. Note that MEGNet and CGCNN are the only models with a positive median. Their hull distance errors are biased towards more frequently predicting thermodynamic instability, explaining why they are closest to getting the overall number of stable structures in the test set right (see cumulative precision/recall in @fig:cumulative-precision-recall).

BOWSR has the largest median error, while Voronoi RF has the largest IQR. Note that MEGNet and CGCNN are the only models with a positive median. Their hull distance errors are biased towards more frequently predicting thermodynamic instability, explaining why they are closest to correctly predicting the overall number of stable structures in the test set (see @fig:cumulative-precision-recall).

### Classification Histograms using Model-Predicted Energies

{#if mounted}
<HistClfPredHullDistModels />
{/if}

> @label:fig:hist-clf-pred-hull-dist-models Distribution of model-predicted hull distance colored by stability classification. Models are sorted from top to bottom by F1 score. The thickness of the red and yellow bands shows how often models misclassify as a function of how far away from the convex hull they place a material. While CHGNet and M3GNet perform almost equally well overall, these plots reveal that they do so via different trade-offs. M3GNet commits fewer false negatives but more false positives predictions compared to CHGNet. In a real discovery campaign, false positives have a higher opportunity cost than false negatives since they result in wasted DFT relaxations or even synthesis time in the lab. A false negative by contrast is just one missed opportunity out of many. This observation is also reflected in the higher TPR and lower TNR of M3GNet vs CHGNet in @fig:metrics-table, as well as the lower rolling MAE for CHGNet vs M3GNet on the stable side (left half) of @fig:rolling-mae-vs-hull-dist-models and vice-versa on the unstable side (right half).

Note the CGCNN+P histogram is more strongly peaked than CGCNN's which agrees better with the actual DFT ground truth [distribution of hull distances](/data#--target-distribution) in our test set. This explains why CGCNN+P performs better as a regressor, but also reveals how it can perform simultaneously worse as a classifier. By moving predictions closer to the stability threshold at 0 eV/atom above the hull, even small errors are significant enough to tip a classification over the threshold.

## Measuring extrapolation performance from WBM batch robustness

As a reminder, the WBM test set was generated in 5 successive batches, each step applying another element replacement to an MP source structure or a new stable crystal generated in one of the previous replacement rounds. The likelihood by which one element replaces another is governed by ISCD-mined chemical similarity scores for each pair of elements. Naively, one would expect model performance to degrade with increasing batch count, as repeated substitutions should on average 'diffuse' deeper into uncharted regions of material space, requiring the model to extrapolate more. We observe this effect for some models much more than others.

@Fig:rolling-mae-vs-hull-dist-wbm-batches-models shows the rolling MAE as a function of distance to the convex hull for each of the 5 WBM rounds of elemental substitution. These plots show a stronger performance decrease for Wrenformer and Voronoi RF than for UIPs like M3GNet, CHGNet, MACE and even force-less GNNs with larger errors like MEGNet and CGCNN.

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

> @label:fig:rolling-mae-vs-hull-dist-wbm-batches-models Rolling MAE as a function of distance to the convex hull for different models. On WBM batch 1, Wrenformer performs comparably to the SOTA UIPs M3GNet and CHGNet. However, the M3GNet and CHGNet UIPs show stronger extrapolative performance, as they barely deteriorate in performance on later batches that move further away from the original MP training distribution. Wrenformer, by contrast, exhibits a pronounced increase in MAE with batch count.

MEGNet and CGCNN both incur a higher rolling MAE than Wrenformer across all 5 batches and across most or all of the hull distance range visible in these plots. However, similar to the UIPs, MEGNet and CGCNN exhibit very little degradation for higher batch counts. The fact that higher errors for later batches are specific to Wrenformer and Voronoi RF suggests that training only on composition and coarse-grained structural features (spacegroup and Wyckoff positions in the case of Wrenformer; coordination numbers, local environment properties, etc. in the case of Voronoi RF) alone is insufficient to learn an extrapolatable map of the PES.

Given its strong performance on batch 1, it is possible that given sufficiently diverse training data, Wrenformer could become similarly accurate to the UIPs across the whole PES landscape at substantially less training and inference cost. However, the loss of predictive forces and stresses may make Wrenformer unattractive for certain applications even then.

### Largest Errors vs. DFT Hull Distance

Given the large variety of models tested, we asked whether any additional insight into the errors can be gained from looking at how the predictions vary between different models.
In @fig:scatter-largest-errors-models-mean-vs-true-hull-dist we see two distinct groupings emerge when looking at the 200 structures with the largest errors.
This clustering is particularly apparent when points are colored by model disagreement.

{#if mounted}
<LargestErrorScatterSelect />
{/if}

> @label:fig:scatter-largest-errors-models-mean-vs-true-hull-dist DFT vs predicted hull distance (average over all models) for the 200 largest error structures colored by model disagreement (as measured by the standard deviation in hull distance predictions from different models) and markers sized by the number of atoms in the structures.
> This plot shows that high-error predictions are biased towards predicting too small hull distances.
> This is unsurprising considering MP training data mainly consists of low-energy structures.
> However, note the clear color separation between the mostly blue low-energy bias predictions and the yellow/red high-error prediction.
> Blue means models are in good agreement, i.e. all models are \"wrong\" together.
> Red/yellow are large-error predictions with little model agreement, i.e. all models are wrong in different ways.
> Some of the blue points with large errors yet good agreement among models may be accurate ML predictions for a DFT relaxation gone wrong.
> Zooming in on the blue points reveals that many of them are large. Larger markers correspond to larger structures where DFT failures are less surprising.
> This suggests ML model committees might be used to cheaply screen large databases for DFT errors in a high-throughput manner.

### MEGNet formation energies from UIP-relaxed structures

{#if mounted}
<MetricsTableMegnetUipCombos />
{/if}

> @label:fig:metrics-table-megnet-uip-combos Except for MEGNet RS2RE, all rows in this table show metrics for the task of predicting the relaxed energy given the unrelaxed structure (IS2RE). Metrics in rows labeled M3GNet→MEGNet and CHGNet→MEGNet are the result of passing M3GNet/CHGNet-relaxed structures into MEGNet for formation energy prediction. Both model combos perform worse than using the respective UIPs on their own with a more pronounced performance drop from CHGNet to CHGNet→MEGNet than M3GNet to M3GNet→MEGnet. This suggests that MEGNet has no additional knowledge of the PES that is not already encoded in the UIPs. However, both combos perform better than MEGNet on its own, demonstrating that ML relaxation imparts some of the latent knowledge of the PES into the relaxed structure from which a more informed final energy prediction is possible. This underscores the utility of ML relaxation at a very low cost for any downstream structure-dependent analysis.

Interestingly, @fig:metrics-table-megnet-uip-combos shows that MEGNet RS2RE performs worse than CHGNet and M3GNet even when it has access to the ground truth DFT relaxed structures when making relaxed energy predictions.
While the regression and classification metrics for MEGNet RS2RE are better than M3GNet→MEGNet and CHGNet→MEGNet, the difference is surprisingly small, especially compared to the regression performance of standalone MEGNet.
This suggests that for downstream use, at least by an imprecise method such as another ML model, the UIP-relaxed structures are almost as helpful as the DFT-relaxed ones.

CHGNet, M3GNet and MEGNet were all trained on slightly different targets.
While MEGNet was trained to predict formation energies, M3GNet and CHGNet were trained on raw DFT energies with the important distinction that CHGNet targets include [MP2020 energy correction scheme](https://pymatgen.org/pymatgen.entries.compatibility.html#pymatgen.entries.compatibility.MaterialsProject2020Compatibility) [@wang_framework_2021] while M3GNet targets do not.

Changing the target from raw DFT energies to formation energies in principle merely constitutes a change of gauge since the difference is just a linear transformation of subtracting the elemental reference energies weighted by composition.
But in practice, one might expect the problem of stability prediction to be easier when trained on formation energies, as the model output is one less step removed from the final target of convex hull distance.
However, empirical evidence so far suggests the opposite.
Both UIP→MEGNet combos perform worse than those same UIPs by themselves.
There are too many confounding effects at play to draw firm conclusions but this tentatively suggests using formation energy as targets offers no benefits over raw DFT.

### Element Prevalence vs Model Error

{#if mounted}
<ElementPrevalenceVsError />
{/if}

> @label:fig:element-prevalence-vs-error The y-axis is the average error of each model on each element averaged over all structures in the test set containing said element weighted by that structure's composition. The x-axis is each element's prevalence in the MP training set, i.e. number of structures containing this element. We observe little correlation between element prevalence and model error. Instead, oxygen, the element with highest representation, is among the higher error elements for most models. Similarly, fluorine, the element with highest average error across models has very good training support at ~12 k samples, suggesting at the sample counts MP provides for its 89 elements, chemistry more than training coverage determines element errors.

### Cumulative MAE

{#if mounted}
<CumulativeMae />
{/if}

> @label:fig:cumulative-mae-rmse CHGNet is consistently the lowest error regressor for materials it predicts as stable. Similar to @fig:cumulative-precision-recall, this figure shows the cumulative MAE over the course of a discovery campaign. As in @fig:cumulative-precision-recall, each model's line starts with those structures the model predicts as most stable and ends with those structures the model predicts as barely stable, i.e. predicted to lie directly on the convex hull.
