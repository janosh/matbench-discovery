<script>
  import { browser } from '$app/environment'
  import { onMount } from 'svelte'
  import { repository as repo } from '$site/package.json'
  import BoxHullDistErrors from '$figs/box-hull-dist-errors.svelte'
  import CumulativeMae from '$figs/cumulative-mae.svelte'
  import CumulativePrecisionRecall from '$figs/cumulative-precision-recall.svelte'
  import EachParityModels from '$figs/each-parity-models-9x2.svelte'
  import ElementPrevalenceVsError from '$figs/element-prevalence-vs-error.svelte'
  import HistClfPredHullDistModels from '$figs/hist-clf-pred-hull-dist-models-9x2.svelte'
  import HullDistParityWrenformerFailures from '$figs/hull-dist-parity-wrenformer-failures.svelte'
  import LargestErrorScatterSelect from './largest-error-scatter-select.svelte'
  import MetricsTableUniqProtos from '$figs/metrics-table-uniq-protos.svelte'
  import MetricsTable from '$figs/metrics-table.svelte'
  import MetricsTableFirst10k from '$figs/metrics-table-first-10k.svelte'
  import MetricsTableMegnetUipCombos from '$figs/metrics-table-uip-megnet-combos.svelte'
  import RocModels from '$figs/roc-models.svelte'
  import RollingMaeVsHullDistModels from '$figs/rolling-mae-vs-hull-dist-models.svelte'
  import RunTimeBars from '$figs/model-run-times-bar.svelte'
  import SpacegroupSunburstWbm from '$figs/spacegroup-sunburst-wbm.svelte'
  import SpacegroupSunburstWrenformerFailures from '$figs/spacegroup-sunburst-wrenformer-failures.svelte'
  import MPvsMPTrjVsWbmArityHist from '$figs/mp-vs-mp-trj-vs-wbm-arity-hist.svelte'
  import RollingMaeVsHullDistWbmBatchesAlignn from '$figs/rolling-mae-vs-hull-dist-wbm-batches-alignn.svelte'
  import RollingMaeVsHullDistWbmBatchesCgcnn from '$figs/rolling-mae-vs-hull-dist-wbm-batches-cgcnn.svelte'
  import RollingMaeVsHullDistWbmBatchesChgnet from '$figs/rolling-mae-vs-hull-dist-wbm-batches-chgnet.svelte'
  import RollingMaeVsHullDistWbmBatchesGnome from '$figs/rolling-mae-vs-hull-dist-wbm-batches-gnome.svelte'
  import RollingMaeVsHullDistWbmBatchesM3gnet from '$figs/rolling-mae-vs-hull-dist-wbm-batches-m3gnet.svelte'
  import RollingMaeVsHullDistWbmBatchesMace from '$figs/rolling-mae-vs-hull-dist-wbm-batches-mace.svelte'
  import RollingMaeVsHullDistWbmBatchesMatterSim from '$figs/rolling-mae-vs-hull-dist-wbm-batches-mattersim.svelte'
  import RollingMaeVsHullDistWbmBatchesOrb from '$figs/rolling-mae-vs-hull-dist-wbm-batches-orb.svelte'
  import RollingMaeVsHullDistWbmBatchesMegnet from '$figs/rolling-mae-vs-hull-dist-wbm-batches-megnet.svelte'
  import RollingMaeVsHullDistWbmBatchesSevenNet from '$figs/rolling-mae-vs-hull-dist-wbm-batches-sevennet.svelte'
  import RollingMaeVsHullDistWbmBatchesVoronoiRf from '$figs/rolling-mae-vs-hull-dist-wbm-batches-voronoi-rf.svelte'
  import RollingMaeVsHullDistWbmBatchesWrenformer from '$figs/rolling-mae-vs-hull-dist-wbm-batches-wrenformer.svelte'
  import wbm_matminer_2d_umap_url from '$figs/wbm-final-struct-matminer-features-2d-umap.png?url'
  import ScatterLargestErrorsModelsMeanVsTrueHullDist from '$figs/scatter-largest-errors-models-mean-vs-true-hull-dist-all.svelte'
  import MpTrjNSitesHist from '$figs/mp-trj-n-sites-hist.svelte'
  import MPTrjEFormHist from '$figs/mp-trj-e-form-hist.svelte'
  import MPTrjForcesHist from '$figs/mp-trj-forces-hist.svelte'
  import MPTrjStressesHist from '$figs/mp-trj-stresses-hist.svelte'
  import MPTrjMagMomsHist from '$figs/mp-trj-magmoms-hist.svelte'
  import HistWbmEForm from '$figs/hist-wbm-e-form.svelte'
  import HistWbmHullDist from '$figs/hist-wbm-hull-dist.svelte'
  import EFormParityModels from '$figs/e-form-parity-models-7x2.svelte'

  let mounted = false
  onMount(() => (mounted = true))
</script>

<summary class="abstract">

Matbench Discovery simulates using machine learning (ML) energy models as pre-filters to DFT in a high-throughput search for stable inorganic crystals.
We address the disconnect between (i) thermodynamic stability and formation energy and (ii) in-domain vs out-of-distribution performance.
Alongside this paper, we publish a Python package to aid with future model submissions and a growing online leaderboard with further insights into trade-offs between various performance metrics.
To answer the question which ML methodology performs best at materials discovery, our initial release explores random forests, graph neural networks (GNN), one-shot predictors, iterative Bayesian optimizers and universal interatomic potentials (UIP).
Ranked best-to-worst by their test set F1 score on thermodynamic stability prediction, we find CHGNet > M3GNet > MACE > ALIGNN > MEGNet > CGCNN > CGCNN+P > Wrenformer > BOWSR > Voronoi fingerprints with random forest.<br>
The top 3 models are all UIPs, which we declare the winning methodology for ML-guided materials discovery, achieving F1 scores of ~0.6 for crystal stability classification and discovery acceleration factors (DAF) of up to 5x on the first 10k most stable predictions compared to dummy selection from our test set.<br>
We also highlight a sharp disconnect between commonly used global regression metrics and more task-relevant classification metrics.
Accurate regressors are susceptible to unexpectedly high false-positive rates if those accurate predictions lie close to the decision boundary at 0 eV/atom above the convex hull where most materials are.
Our results highlight the need to for the field to focus more on classification than regression metrics, since the former actually correlate with improved stability hit rate.
Finally, we share valuable insights for maintainers of high throughput materials databases by demonstrating that universal potentials have matured enough to play a useful role as triaging tools for effectively allocating compute budget in high-throughput DFT.

</summary>

## Introduction

Material science can be viewed as a combinatorial problem of mixing and arranging different atoms to leverage the complex range of properties that emerge.
To date, $~10^5$ combinations have been tested experimentally
[@bergerhoff_inorganic_1983; @belsky_new_2002]
and $~10^7$ have been simulated
[@jain_commentary_2013; @saal_materials_2013; @curtarolo_aflow_2012; @draxl_nomad_2018].
[@davies_computational_2016] identified
$~10^{10}$ possible quaternary materials allowed by electronegativity and charge-balancing rules.
The space of quinternaries and higher is even less explored, leaving vast numbers of potentially useful materials to be discovered.
The discovery of new materials is a key driver of technological progress and lies on the path to more efficient solar cells, lighter and longer-lived batteries, smaller and more efficient transistor gates just to name a few.
In light of our sustainability goals, these advances cannot come fast enough. Any speed-up new discovery methods might yield should be leveraged to their fullest extent.

Despite significant advances in empirical, theoretical and computational materials science, discovering new materials still requires complex calculations, labor-intensive trial-and-error experimentation, and often happens fortuitously rather than through rational design.
Machine learning (ML) methods efficiently extract and distill trends from huge datasets, and can handle high dimensionality, multiple objectives [@riebesell_pushing_2024], uncertainty [@borg_quantifying_2022; @goodall_rapid_2022; @zhu_fast_2023], and noisy or sparse data [@depeweg_decomposition_2017; @bartel_critical_2020], making them powerful additions to the computational materials science tool set.

ML models are less accurate and reliable but orders of magnitude faster than ab initio simulation.
This makes them most suitable for use in high-throughput (HT) searches to triage more expensive, higher-fidelity simulation methods.
The use of neural networks for learning the Kohn-Sham density-functional theory (DFT) potential energy surface (PES) can be traced as far back as [@behler_generalized_2007].
This work kicked off rapid advances and significant efforts to fit ever more sophisticated ML models to known samples of the PES.
Initially, most of these models were trained and deployed as interatomic potentials (also known as force fields) to study known materials of interest, a workflow that requires curating custom training data for each new system of interest [@bartok_machine_2018; @deringer_general-purpose_2020].
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
   Distance to the convex hull spanned by competing phases in the same chemical system is a much more informative indicator of thermodynamic stability.
   Moreover, ML models relying on relaxed crystal structures as input render any discovery pipeline circular since obtaining relaxed structures requires computationally expensive DFT simulations, thereby depending on the very process we intend to accelerate.

2. **Opportunity cost**: Accurate regressors are susceptible to unexpectedly high false-positive rates if those accurate predictions lie close to the decision boundary.
   Looking purely at global metrics like $\text{MAE}$, $\text{RMSE}$ and $R^2$ can give practitioners a false sense of security about their model's reliability.
   Failed experiments incur a high opportunity cost by wasting lab time and resources.

3. **Scalability**: Future discovery efforts are likely to encounter large data regimes.
   Small benchmarks can lack chemical diversity, obfuscate poor scaling relations or poor out-of-distribution (OOD) performance (see @fig:wbm-final-struct-matminer-features-2d-umap for a UMAP projection of train and test sets).
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
Our findings show that UIPs outperform all other methodologies we tested both in accuracy and robustness.

### Evaluation Framework for Materials Discovery

We propose a novel evaluation framework that places no constraints on the type of data a model trains on as long as it would be available to a practitioner conducting a real materials discovery campaign.
We only enforce that at test time, all models must make predictions on the convex hull distance of the relaxed structure with only the unrelaxed structure as input.
This setup avoids circularity in the discovery process, as unrelaxed structures can be cheaply enumerated through elemental substitution methodologies and do not contain information inaccessible in a prospective discovery campaign.
We choose to focus on the relaxed structure's convex hull distance as a measure of thermodynamic stability rather than the formation energy as it informs the decision on whether to pursue a potential candidate crystal.
This decision was also motivated by [@bartel_critical_2020] who found even composition-only models capable of predicting DFT formation energies with useful accuracy.
However, when tasking those same models with predicting decomposition enthalpy, performance deteriorated sharply.
This insight meant that ML models are much less useful than DFT for discovering new inorganic solids.
Moreover, the qualitative leap in performance from Roost [@goodall_predicting_2020], the best compositional model benchmarked, to CGCNN [@xie_crystal_2018], the single structural model they tested, highlights that structure plays a crucial role in determining material stability.
Hence our decision to restrict test-time model input to unrelaxed structures to avoid the above-mentioned circularity and measure true prospective utility while still allowing for subtle atomic-configuration-informed energy estimates.
We recognize that our ground truth is ignorant of entropic stabilization and metastable states.
While these factors influence experimental synthesizability, we have to neglect them when simulating an HT search since the convex hull distance as predicted by zero-temperature DFT is the best proxy of true crystal stability available for large datasets.

Standard practice in ML benchmarks is to hold all variables fixed -- most importantly the training data -- and vary only the model architecture to isolate architectural effects on the performance.
We deliberately deviate from this practice due to diverging objectives from common ML benchmarks.
Our goal is to identify the best methodology for accelerating materials discovery.
What kind of training data a model can ingest is part of its methodology. Unlike energy-only models, UIPs benefit from the additional training data provided by the forces and stresses recorded in DFT relaxations. This allows them to learn a fundamentally higher fidelity model of the physical interactions between ions. That is a genuine advantage of the architecture and something any benchmark aiming to identify the optimal methodology for materials discovery must reflect, as ours does.
As a result, our benchmark contains models trained on varying data sets.

We define the Materials Project (MP) [@jain_commentary_2013] database release ([v2022.10.28](https://docs.materialsproject.org/changes/database-versions) as the maximum allowed training set for any valid model submission. Any subset of the energies, forces, stresses, magnetic moments or any other recorded DFT data are valid training targets.
All of these would be available to a practitioner performing a real materials discovery campaign and hence are permitted for any model submission.
Models may train on the complete set of relaxation frames, or any subset thereof such as the final relaxed structures only. Our test set consists of the unrelaxed structures in the WBM dataset [@wang_predicting_2021]. Their target values are the PBE formation energies of the corresponding DFT-relaxed structures.

#### Materials Project Training Set

The Materials Project is a well-known effort to calculate the properties of all inorganic materials using high-throughput ab initio methods.
Seeded from a subset of the Inorganic Crystal Structure Database (ICSD) [@allen_crystallographic_1999], the initial release of the database consisted of ~9 k crystals.
At the time of writing, the Materials Project database [@jain_commentary_2013] has grown to ~154 k crystals, covering diverse chemistries and providing relaxed and initial structures as well as the relaxation trajectory for every entry.

Our benchmark defines the training set as all data available from the [v2022.10.28 MP release](https://docs.materialsproject.org/changes/database-versions#v2022.10.28).
We recorded a snapshot of energies, forces, stresses and magnetic moments for all MP ionic steps on 2023-03-15 as the canonical training set for Matbench Discovery, and provide convenience functions through our [Python package](https://pypi.org/project/matbench-discovery) for easily feeding that data into future model submissions to our benchmark.

The flexibility in specifying the maximum allowed data set is intended to allow authors to experiment with and exploit the large variety of available data.
This choice is motivated by two factors.
First, whilst it may appear that models that train on multiple snapshots containing energies, forces and stresses are getting more training data than models that are only trained or able to be trained on the energies of the relaxed structures the critical factor is that all this additional data was generated as a byproduct of the workflow to produce the relaxed structures i.e. all models have been trained from the same total cost of training data acquisition.
If some architectures or approaches can leverage more of this byproduct data to make improved predictions this is a fair comparison between the two models.
This approach diverges philosophically from other benchmarks such as the OCP and Matbench where it has been more common to subcategorize different models and look at multiple alternative tasks (e.g. composition-only vs structure-available in Matbench or IS2RS, IS2RE, S2EF in OCP) and do not make direct comparisons of this manner.
Second, recent work in the space from [@li_critical_2023; @li_exploiting_2023] has claimed that much of the data in large databases like MP are redundant and that models can be trained more effectively by taking a subset of these large data pools.
We believe that from a systems-level perspective identification of novel cleaning and or active-learning strategies that do not leak signal from the test set are equally as important as architectural improvements as they can lead to similar improvements in performance, particularly given the prevalence of errors in high-throughput DFT that can disrupt learning of the PES.
Consequently, such strategies where they lead to improved performance should be able to be recognized within the benchmark.
We encourage submissions to submit ablation studies showing how different system-level choices affect performance.

We highlight for example data sets that are valid within the rules of the benchmark that take advantage of these freedoms.
The first is the MPF.2021.2.8 data set [@chen_universal_2022] curated to train the M3GNet model which takes a subset of materials from the v2021.02.08 MP release.
Rather than taking every ionic step from the relaxation trajectory this data set opts to select only the initial, final and one intermediate structure for each material to avoid biasing the data set towards examples where more ionic steps were needed to relax the structure.
Additionally, the curators of the MPF.2021.2.8 data set down-sampled the v2021.02.08 release significantly to select a subset of calculations that they believed to be most self-consistent.
The MPF.2021.2.8 is a proper subset of the training data as no materials were deprecated between the v2021.02.08 and v2022.10.28 database releases.
The second data set we highlight, with which several of the UIP models have been trained, is the MPtrj dataset [@deng_chgnet_2023]. This data set was curated from the earlier v2021.11.10 MP release. The MPtrj dataset is a proper subset of the allowed training data but several potentially anomalous examples from within MP were cleaned out of the data set before the frames were subsampled to remove. It is worth noting that the v2022.10.28 release contains a small number of additional Perovskite structures not found in MPtrj that could be added to the training set within the scope of the benchmark.

We note that the v2023.11.1 deprecated a large number of calculations so data queried from subsequent database releases is not considered valid for this benchmark.

#### WBM Test Set

The WBM data set [@wang_predicting_2021] consists of ~257 k structures generated via chemical similarity-based elemental substitution of MP source structures followed by DFT relaxation and calculating each crystal's convex hull distance.
Which element replaces an existing one in a given source structure was determined by random sampling according to the weights in a chemical similarity matrix data-mined from the ICSD [@glawe_optimal_2016].

The WBM authors performed 5 iterations of this substitution process (we refer to these as batches).
After each step, the newly generated structures found to be thermodynamically stable after DFT relaxation flow back into the source pool to partake in the next round of substitution.
This split of the data into batches of increasing substitution count is a unique and compelling feature of the test set as it allows out-of-distribution (OOD) testing by seeing if model performance degrades with substitution count.
A higher number of elemental substitutions on average carries the structure further away from the region of material space covered by the MP training set (see @fig:rolling-mae-vs-hull-dist-wbm-batches-models, @fig:wbm-final-struct-matminer-features-2d-umap for details).
Whilst this batch information makes the WBM data set an exceptionally useful data source for examining the extrapolation performance of ML models, our evaluation primarily looks at metrics that consider all batches as a single test set.

Throughout this work, we define stability as being on or below the convex hull of the MP training set ($E_\text{MP hull dist} \leq 0$).
~42k out of ~257k materials in WBM satisfy this criterion.
Of these, ~33k are unique prototypes, meaning they have no matching structure prototype in MP nor another higher-energy duplicate prototype in WBM.
Our code treats the stability threshold as a dynamic parameter, allowing for future model comparisons at different thresholds.
This may reveal systematic differences in stability calibration, as hinted by the different TPR/TNR tradeoffs observed for CHGNet and M3GNet in @tab:metrics-table-first-10k.
A looser stability threshold (e.g. $E_\text{MP hull dist} \leq 0.05 eV/atom$) might benefit a model like M3GNet which is more conservative in its stability predictions, whereas a stricter threshold may benefit a model like CHGNet that is more optimistic in its stability predictions.
For initial analysis in this direction, see @fig:roc-models in the SI.

As WBM explores regions of materials space not well sampled by MP, many of the discovered materials that lie below MP's convex hull are not stable relative to each other.
Of the ~42k that lie below the MP convex hull less than half, or around ~20k, remain on the joint MP+WBM convex hull.
This observation suggests that many WBM structures are repeated samples in the same new chemical spaces.
It also highlights a critical aspect of this benchmark in that we purposely operate on an incomplete convex hull.
Only current knowledge of competing points on the PES is accessible to a real discovery campaign and our metrics are designed to reflect this.

### Models

To test a wide variety of methodologies proposed for learning the potential energy landscape, our initial benchmark release includes 10 models.

1. **MACE** [@batatia_mace_2023] (UIP-GNN) - MACE builds upon the recent advances [@batzner_equivariant_2022; @thomas_tensor_2018] in equivariant neural network architectures by proposing an approach to computing high-body-order features efficiently via Atomic Cluster Expansion [@drautz_atomic_2019].
   Unlike the other UIP models considered MACE was primarily developed for molecular dynamics of single material systems and not the universal use case studied here.
   It is the only $E(3)$-equivariant model we tested.

2. **CHGNet** [@deng_chgnet_2023] (UIP-GNN) - CHGNet is a UIP for charge-informed atomistic modeling.
   Its distinguishing feature is that it was trained to predict magnetic moments on top of energies, forces and stresses in the MPtrj dataset consisting of relaxation trajectories for ~1.5 million MP structures (see @sec:eda for detailed analysis of this dataset).
   By modeling magnetic moments, CHGNet learns to accurately represent the orbital occupancy of electrons which allows it to predict both atomic and electronic degrees of freedom.

3. **M3GNet** [@chen_universal_2022] (UIP-GNN) - M3GNet is a GNN-based universal interatomic potential (UIP) (as in covering the full periodic table) for materials trained on up to 3-body interactions in the initial, middle and final frame of MP DFT relaxations.
   The model takes the unrelaxed input and emulates structure relaxation before predicting energy for the pseudo-relaxed structure.

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
<MetricsTableUniqProtos />
{/if}

> @label:tab:metrics-table-uniq-protos Classification and regression metrics for all models tested on our benchmark ranked by F1 score.
> The heat map ranges from yellow (best) to blue (worst) performance.
> DAF = discovery acceleration factor is the ratio of model precision to percentage of stable structures in the test set, TPR = true positive rate = $\frac{\text{TP}}{\text{TP} + \text{FN}}$, TNR = true negative rate = $\frac{\text{TN}}{\text{TN} + \text{FP}}$, MAE = mean absolute error, RMSE = root mean squared error.
> The dummy classifier uses the `scikit-learn` `stratified` strategy of randomly assigning stable/unstable labels according to the training set prevalence.
> The dummy regression metrics MAE, RMSE and $R^2$ are attained by always predicting the test set mean.
> The three universal interatomic potentials MACE, CHGNet and M3GNet take the top 3 spots by F1 score and achieve an impressive performance gap esp. in regression metrics to the other 7 force-free models, highlighting the high information content of forces and stresses for learning a high-fidelity PES.
> Voronoi RF, CGCNN and MEGNet perform worse than dummy in regression metrics but better than dummy on some classification metrics which are more indicative of model utility for guiding discovery, demonstrating that regression metrics alone can be misleading.
> This increases 'OOD-ness' of WBM, this table is restricted to the subset of materials with no matching structure prototype in MP nor another higher-energy duplicate prototype in WBM, reducing the test set size to $\text{257}\to\text{215.5k}$.
> See SI @tab:metrics-table for metrics on the whole WBM test set and details on prototype matching.

@tab:metrics-table-uniq-protos shows performance metrics for all models included in the initial release of Matbench Discovery.
MACE takes the top spot and emerges as the current SOTA for ML-guided materials discovery.
It outperforms the 9 other models on all 6 classification metrics but is slightly overtaken by CHGNet in 2 of 3 regression metrics, $R^2$ and RMSE, which are sensitive to outliers.
The discovery acceleration factor (DAF) measures how many more stable structures a model found compared to the dummy discovery rate achieved by randomly selecting test set crystals, formally the DAF is the ratio of the precision to the prevalence.
The maximum possible DAF is the inverse of the prevalence or dummy discovery rate which on our dataset is $(\text{33k} / \text{215k})^{-1} \approx 6.5$.
Thus the current SOTA of 3.85 achieved by MACE leaves room for improvement.
However, evaluating each model on the much smaller subset of its 10k most stable predictions as shown in @tab:metrics-table-first-10k, CHGNet takes the lead in all 9 metrics, in particular achieving an impressive DAF of 5.29 which approaches optimal performance.

We find a large performance gap between models that make one-shot predictions directly from unrelaxed inputs (MEGNet, Wrenformer, CGCNN, CGCNN+P, ALIGNN, Voronoi RF) compared to UIPs that learn from force and stress labels to emulate DFT relaxation (CHGNet, M3GNet, MACE).
While the F1 scores and DAFs of force-free models are surprisingly less affected, their regression metrics such as $R^2$and RMSE are significantly worse.
Of the force-free models, only ALIGNN, BOWSR and CGCNN+P achieve positive predictive power $R^2$.
Negative $R^2$ means model predictions explain the observed variation in the data less than simply predicting the test set mean.
In other words, these models are not predictive in a global sense (across the full dataset range).
However, even models with negative $R^2$ can be locally good in the positive and negative tails of the test set distribution.
They suffer most in the mode of the distribution near the stability threshold of 0 eV/atom above the hull.
This reveals an important shortcoming of $R^2$ as a metric for classification tasks like stability prediction.

The results for M3GNet depart from the trend that F1 is broadly rank-correlated with the other classification metrics.
M3GNet achieves a high true positive rate ($\text{TPR} = 0.8$) but a notably lower true negative rate ($\text{TNR} = 0.81$) than CHGNet and ALIGNN with the closest F1 scores.
@fig:rolling-mae-vs-hull-dist-models provides a visual explanation for this observation.
M3GNet has a lower rolling mean of the absolute errors (rolling MAE) as a function of hull distance for materials above the convex hull (see right half of plot) but incurs comparably large errors for materials below the hull (left half of plot).
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

> @label:fig:cumulative-precision-recall This figure measures model utility for materials discovery campaigns of varying sizes by plotting the precision and recall as a function of the number of model predictions validated.
> CHGNet initially achieves the highest cumulative precision and recall.
> However, for campaigns that validate ~25k+ model predictions, MACE overtakes CHGNet in both precision and recall and maintains that lead till the end.
> A typical discovery campaign will rank hypothetical materials by model-predicted hull distance from most to least stable and validate the most stable predictions first.
> A higher fraction of correct stable predictions corresponds to higher precision and fewer stable materials overlooked correspond to higher recall.
> This figure highlights how different models perform better or worse depending on the length of the discovery campaign.
> The UIPs (CHGNet, M3GNet, MACE) are seen to offer significantly improved precision on shorter campaigns of ~20k or less materials validated as they are less prone to false positive predictions among highly stable materials.

@Fig:cumulative-precision-recall has models rank materials by model-predicted hull distance from most to least stable; materials furthest below the known hull at the top, materials right on the hull at the bottom.
For each model, we iterate through that list and calculate at each step the precision and recall of correctly identified stable materials.
This simulates exactly how these models would be used in a prospective materials discovery campaign and reveals how a model's performance changes as a function of the discovery campaign length. As a practitioner, you have a certain amount of resources available to validate model predictions. These curves allow you to read off the best model given these conditions and based on the optimal trade-off between fewer false positives (precision) or fewer negatives (recall) for the discovery task at hand.
In this case, it so happens that CHGNet achieves the highest precision _and_ recall at any number of screened materials.

A line terminates when a model believes there are no more materials in the WBM test set below the MP convex hull.
The dashed vertical line shows the actual number of stable structures in our test set.
All models are biased towards stability to some degree as they all overestimate this number, most of all BOWSR by 133%, least of all MEGNet by 30%.
This is only a problem in practice for exhaustive discovery campaigns that validate _all_ stable predictions from a model.
More frequently, model predictions will be ranked most-to-least stable and validation stops after some pre-determined compute budget is spent, say, 10k DFT relaxations.
In that case, the concentration of false positive predictions that naturally accumulates near the less stable end of the candidate list can be avoided with no harm to the campaign's overall discovery rate (see @tab:metrics-table-first-10k a similar metrics table considering only the 10k materials predicted by each model to be furthest below the known convex hull).

The diagonal Optimal Recall line would be achieved if a model never made a false negative prediction and stopped predicting stable crystals exactly when the true number of stable materials is reached.
Zooming in on the top-left corner of the precision plot, we observe that CHGNet is the only model without a sudden drop in precision right at the start of the discovery campaign.
This means CHGNet is the only model whose first few hundred most stable materials are largely actually stable.
MACE in particular suffers from several initial failure cases.
MACE underpredicts the energy of these unstable materials by several eV/atom (see @fig:each-parity-models), resulting in MACE's initial precision starting at 0 before abruptly climbing to ~0.8.
We experienced this behavior with multiple independently trained MACE variants; even more so with a checkpoint we received directly from the MACE authors which was trained on the M3GNet training set. The results shown here are from a superior MACE we trained ourselves on the much larger MPtrj dataset.
Even M3GNet, the other universal potential and 2nd best model, exhibits the early-on precision drop.
As a result, CHGNet has a strong lead over M3GNet until reaching ~3k screened materials.
From there, CHGNet and M3GNet slowly converge until they almost tie at a precision of 0.52 after ~56k screened materials.
At that point, CHGNet's list of stable predictions is exhausted while M3GNet continues, dropping in precision to 0.45 at 76 k, attributable to many false positives near the end of the list of stable predictions.

All force-free models exhibit much stronger early-on precision drop, falling to 0.6 or less in the first 5k screened materials. Many of these models (all except BOWSR, Wrenformer and Voronoi RF) display an interesting hook shape in their cumulative precision, recovering again slightly in the middle of the simulated campaign between ~5k and up to ~30k before dropping again until the end.

{#if mounted}
<RollingMaeVsHullDistModels />
{/if}

> @label:fig:rolling-mae-vs-hull-dist-models Universal potentials are more reliable classifiers because they exit the red triangle earliest.
> These lines show the rolling MAE on the WBM test set as the energy to the convex hull of the MP training set is varied.
> Lower is better.
> Inside the large red 'triangle of peril', models are most likely to misclassify structures.
> As long as a model's rolling MAE remains inside the triangle, its mean error is larger than the distance to the convex hull.
> If the model's error for a given prediction happens to point towards the stability threshold at $E_\text{above MP hull} = 0$, its average error will change the stability classification from true positive/negative to false negative/positive.
> The width of the 'rolling window' box indicates the width over which prediction errors were averaged.

@Fig:rolling-mae-vs-hull-dist-models provides a visual representation of the reliability of different models as a function of a material's DFT distance to the MP convex hull.
The lines show the rolling mean absolute error (RMAE) of model-predicted hull distances vs DFT.
The red-shaded area, which we coin the 'triangle of peril', emphasizes the zone where the average model error surpasses the distance to the stability threshold at 0 eV.
As long as the rolling MAE remains within this triangle, the model is most susceptible to misclassifying structures.
Because the average error is larger than the distance to the classification threshold at 0, it is large enough to flip a correct classification into an incorrect one (if the error happens to point toward the stability threshold).

The sooner a model exits the triangle on the left side (negative DFT hull distance), the less likely it is to incorrectly predict stable structures as unstable, thereby reducing false negatives.
Exciting early on the right side (positive DFT hull distance) results in a lower likelihood of predicting unstable structures as stable, leading to fewer false positives.

On both sides, MACE exits the triangle first at a remarkably low $~\pm30 meV/atom$ hull distance followed closely by CHGNet, consistent with the high accuracy of these models in @tab:metrics-table-uniq-protos.
M3GNet performs similarly to CHGNet for stable structures but has a widening performance gap for increasingly unstable structures on the right side of the plot.
This leads to M3GNet's relatively low TNR of 0.81 compared to CHGNet's 0.87 in @tab:metrics-table-uniq-protos.

All models tend to have lower rolling error towards the plot's left edge compared to the right edge.
This imbalance shows that models are more prone to false positive than false negative predictions. In other words, all models are less likely to predict a material at −0.2 eV/atom DFT hull distance as unstable than a material at 0.2 eV/atom DFT hull distance as stable.
From a practitioner's standpoint, this is undesirable as attempting to validate an unstable material usually has a much higher opportunity cost than missing a stable one.
We hypothesize this error imbalance is due to the MP training set having an uncharacteristically high share of stable materials, causing statical models trained on it to be biased towards low-energy predictions even for very energy atomic configurations.

As the convex hull becomes more thoroughly sampled by future discovery, the fraction of unknown stable structures decreases, naturally leading to less enriched future test sets.
This has several implications, including allowing for higher maximum DAFs but also skewing the hull distance distribution towards positive values (i.e. the right half of @fig:rolling-mae-vs-hull-dist-models) as there are fewer undersampled chemical spaces. Consequently, to have accurate stability classifications, model accuracy will need to improve concomitantly.
For @fig:rolling-mae-vs-hull-dist-models, this means models will need to become more accurate to still be able to exit the shaded triangle on the left side of the plot.

## Discussion

We have demonstrated the effectiveness of ML-based triage in HT materials discovery and posit that the benefits of including ML in discovery workflows now clearly outweigh the costs.
@tab:metrics-table-uniq-protos shows in a realistic benchmark scenario that several models achieve a discovery acceleration greater than 2.5 across the whole dataset and up to 5 when considering only the 10k most stable predictions from each model (@tab:metrics-table-first-10k).
When starting this project, we were unsure which is the most promising ML methodology for HT discovery.
Our findings demonstrate a clear superiority in the accuracy and extrapolation performance of UIPs like CHGNet, M3GNet and MACE.
Modeling forces enables these models to chart a path through atomic configuration space closer to the DFT-relaxed structure from where a more informed final energy prediction is possible.
The consistently linear log-log learning curves observed in related literature [@vonlilienfeld_retrospective_2020] suggest that further decreases in the error of UIPs can be readily unlocked with increased training data.

Despite its fast pace of progress, the field of universal ML potentials is still in its early stages.
To realize the full potential of these models, we will need to pivot significant resources to generate large quantities of higher-than-PBE fidelity training data.
That said, we believe our work (e.g. @fig:rolling-mae-vs-hull-dist-wbm-batches-models) and others [@lan_adsorbml_2023] have caught early glimpses that sufficiently large training sets allow us to narrow the performance gap between in-domain and out-of-domain performance.
This should further increase the discovery acceleration factor, resulting in a virtuous cycle between ML models accelerating the creation of new training data while becoming more accurate and even stronger accelerators after retraining.
Meanwhile, the setup cost and learning curve of ML potentials are steadily decreasing, likely with a much lower floor than DFT given the intricacies of HT ab-initio simulation.
Combined with their already strong performance, this suggests that UIPs combined with data accumulation over the next few years may unlock our way to a comprehensive yet cheap map of the PES that truly commoditizes HT materials discovery.

Despite this optimistic outlook, not only has the predictive power UIPs attained for crystal stability over the last 3 years since [@bartel_critical_2020] sobering 2020 analysis yet to enter common knowledge, but many remaining shortcomings invite further research.
For example, all models are biased towards stability, overestimating the correct number of stable materials in our test set anywhere between 30% and 133%.
Moreover, even the best models in our analysis still exhibit a significant fraction of false positive predictions, even for materials as high as 50 meV/atom above the hull.
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

We thank Rickard Armiento, Felix A. Faber and Abhijith S. Parackal for helping develop the evaluation procedures for Wren upon which this work builds. We also thank Rokas Elijosius for assisting in the initial implementation of Wrenformer.

We would like to thank Jason Blake Gibson, Shyue Ping Ong, Chi Chen, Tian Xie, Peichen Zhong and Ekin Dogus Cubuk for helpful discussions.

## Author Contributions

Janosh Riebesell: Methodology, software, data curation, training and testing models, formal analysis. Rhys Goodall: Conceptualization, software, formal analysis. Philipp Benner: Software, training ALIGNN + MACE. Yuan Chiang: Training MACE, data analysis. Bowen Deng: Curating MPtrj dataset and training CHGNet. Alpha Lee: Supervision. Anubhav Jain: Supervision. Kristin Persson: Funding Acquisition.

## Code availability

We welcome further model submissions to our GitHub repository <https://github.com/janosh/matbench-discovery>.

## Data availability

We chose the latest Materials Project (MP) [@jain_commentary_2013] database release ([v2022.10.28](https://docs.materialsproject.org/changes/database-versions) at time of writing) as the training set and the WBM dataset [@wang_predicting_2021] available at <https://figshare.com/articles/dataset/22715158> as the test set for this benchmark.
A snapshot of every ionic step including energies, forces, stresses and magnetic moments in the MP database is available at <https://figshare.com/articles/dataset/23713842>.

## Supplementary Information

### Metrics on full test set and for 10k materials predicted most stable

Unlike @tab:metrics-table-uniq-protos which evaluates model performance on the subset of unique WBM structure prototypes, @tab:metrics-table includes the full WBM test set of 257k materials.
The 44k additional materials in @tab:metrics-table excluded in @tab:metrics-table-uniq-protos comprise 11175 discarded due to having a matching structure prototype in MP plus another 32784 materials which are prototype duplicates of another WBM material with lower energy.
The most noteworthy difference between the two tables is a drop in DAF for all models in @tab:metrics-table compared to @tab:metrics-table-uniq-protos.
MACE for example achieves a DAF of 3.5 on the full test set (@tab:metrics-table) compared to 3.85 on the subset of 215.5k materials with unique prototypes (@tab:metrics-table-uniq-protos).
This ~10% increase is largely due to a ~10% decrease in the fraction of materials below the MP convex hull: 15.3% (32,942 out of 215,488) in @tab:metrics-table-uniq-protos vs 16.7% (42,825 out of 256,963) in the full dataset (@tab:metrics-table).
Since DAF is the ratio of the model's precision for stability prediction to the prevalence of stable structures, a lower prevalence results in a higher DAF.

While prototypes are non-trivial to match, thus potentially introducing bias into the test set by trying to deduplicate them, we still opted to feature @tab:metrics-table-uniq-protos in the main text since the removal of overlapping prototypes with MP makes it more closely reflect a model's true ability to extrapolate to out-of-domain materials.
Prototypes were matched based on an Aflow-style [@hicks_aflow_2021] Wyckoff representation as implemented in `aviary.wren.utils.get_protostructure_label_from_spglib` [@goodall_rapid_2022]. This string encodes the crystal's prototype and is invariant to relaxation (unless the system changes spacegroup during relaxation, usually by gaining symmetry) and will count structures with different lattice parameters as duplicates if both are expected to relax to the same ground state.

{#if mounted}
<MetricsTable />
{/if}

> @label:tab:metrics-table Same as @tab:metrics-table-uniq-protos but computed for all 257k structures in the WBM test set, no duplicate prototypes excluded.
> The prevalence of stable structures is 16.7% (42,825 out of 256,963), placing a ceiling of 6 on the maximum possible DAF.

A real-world discovery campaign is unlikely to validate all stable predictions from a given model as we did in @tab:metrics-table-uniq-protos.
Presumably, it will rank model predictions from most to least stable and follow that list as far as time and compute budget permits.
Assuming that increases in compute resources will allow average future discovery campaigns to grow in scope, we believe 10k model validations to be a reasonable scope for average campaigns.
This is what @tab:metrics-table-first-10k simulates by calculating classification and regression metrics for the 10k test set materials predicted to be most stable by each model.
We again show dummy performance in the bottom row. Note that each model is now evaluated on a different slice of the data. However, the bottom row still shows dummy performance across the whole dataset.
CHGNet and M3GNet achieve an impressive precision of 0.83 and 0.8, respectively.
In concrete terms, this means in a discovery campaign that validates 10 k model predictions from a search pool of 257 k crystals that are chemically dissimilar from the training set and of which 0.167 are stable, CHGNet and M3GNet would deliver 4 stable structures for every 5 predictions validated.
We emphasize that in light of the significant resulting increase in stability hit rate, these models are well worth integrating into future materials searches.

{#if mounted}
<MetricsTableFirst10k />
{/if}

> @label:tab:metrics-table-first-10k Stability prediction metrics for the 10k materials predicted to be most stable by each model.
> Every model picks out a different slice of the test set as the 10k materials it believes to be furthest below the known convex hull.
> CHGNet achieves an impressive F1 score of 0.91 and a discovery acceleration factor (DAF) of 5.29, both approaching the optimal values of 1 and 6, respectively.
> CHGNet notably outperforms MACE on its first 10k predictions sorted by stability which are the most important ones when guiding a discovery effort as these are usually the first materials to be validated by higher fidelity methods.
> This is reflected in the left panel of @fig:cumulative-precision-recall and the top row of @fig:each-parity-models where MACE, unlike CHGNet, suffers a sudden drop in precision right at the start of the discovery campaign from severe underpredictions of a small number of compounds.
> These can be seen as points along the negative $y$ axis in the MACE parity plot in @fig:each-parity-models which are not present for CHGNet.

### ROC Curves

A material is classified as stable if the predicted $E_\text{above hull}$ is below the stability threshold. Since all models predict $E_\text{form}$ (and M3GNet predicted energies are converted to formation energy before stability classification), they are insensitive to changes in the threshold $t$. The receiver operating characteristic (ROC) curve for each model is plotted in @fig:roc-models. M3GNet wins in the area under curve (AUC) with 0.87, coming in 1.34x higher than the worst model, the Voronoi Random Forest. The diagonal 'No skill' line shows the performance of a dummy model that randomly ranks material stability.

{#if mounted}
<RocModels />
{/if}

> @label:fig:roc-models Receiver operating characteristic (ROC) curve for each model. TPR/FPR = true/false positive rate. FPR on the $x$ axis is the fraction of unstable structures classified as stable. TPR on the $y$ axis is the fraction of stable structures classified as stable.

### Parity Plots

{#if mounted}
<EachParityModels />
{/if}

> @label:fig:each-parity-models Parity plots of model-predicted energy distance to the convex hull (based on their formation energy predictions) vs DFT ground truth, color-coded by log density of points.
> Models are sorted left to right and top to bottom by MAE.
> For parity plots of formation energy predictions, see @fig:e-form-parity-models.

@Fig:each-parity-models shows that all models do well for materials far below the convex hull (left side of the plot). Performance for materials far above the convex hull is more varied with occasional underpredictions of the energy of materials far above the convex hull (right side). All models suffer most in the mode of the distribution at $x = 0$.

Two models stand out as anomalous to the general trends.

Wrenformer is the only model with a large number of severe energy overpredictions at $x = 0$ along the positive $y$ axis. We investigate these failure cases in more detail in @fig:wrenformer-failures and find these overpredictions to be dominated by spacegroup 71 with poor representation in the training data.

The other anomalous model is MACE with several severe underpredictions at $x = 0$ along the negative $y$ axis. We investigated these points for common traits in composition or crystal symmetry but noticed no pattern.

Beyond these MACE outliers visible in the plot, MACE exhibits another rare but reproducible type of failure case, in which the final predicted energy after relaxation is off by several orders of magnitude. The largest 'derailed' prediction was $-10^{22}$ eV/atom for `wbm-3-31970` (formula H$_2$Ir). In each case, the MACE relaxation exhausted the maximum number of ionic steps set to 500 and caused a volume implosion from initial cell volumes of hundreds to relaxed cell volumes of tens of cubic Angstrom. Using the checkpoint trained on the M3GNet dataset which we received from the MACE authors, this failure mode occurred for several hundred of the 250k test set crystals. Using the checkpoint we trained ourselves on the MPtrj dataset, it affects only 44 test crystals, suggesting that these holes in the MACE PES can perhaps be fully plugged by further increasing the training set or even changing the loss function.
Further analysis is ongoing.
Since these derailed values are easily identified in practice when actually performing a prospective discovery campaign, we excluded them from the MACE parity plat and all other downstream analyses.

### Spacegroup prevalence in Wrenformer failure cases

{#if mounted}

  <div style="display: flex; gap: 1em; justify-content: space-around; flex-wrap: wrap; margin-bottom: 2em;">
    <SpacegroupSunburstWrenformerFailures />
    <SpacegroupSunburstWbm />
  </div>
  <HullDistParityWrenformerFailures />
{/if}

> @label:fig:wrenformer-failures Symmetry analysis of the 941 Wrenformer failure cases in the shaded rectangle defined by $E_\text{DFT hull dist} < 1$ and $E_\text{ML hull dist} > 1$. Sunburst plot of spacegroups shows that close to 80% of severe energy overestimations are orthorhombic with spacegroup 71. The table on the right shows occurrence counts of exact structure prototypes for each material in the sunburst plot as well as their corresponding prevalence in the training set.

@Fig:wrenformer-failures shows 456 + 194 (~70%) of the failure cases in the shaded rectangle are two prototypes in spacegroup 71.
The occurrence of those same prototypes in the MP training set shows almost no data support for the failing prototypes.
This suggests the reason Wrenformer fails so spectacularly on these structures is that it cannot deal with structure prototypes it has not seen at least several hundred examples of in its training data.
This suggests that there are stronger limitations on how much the discrete Wyckoff-based representation can extrapolate to new prototypes compared to the smooth local-environment-based inputs to GNN-type models.

### Hull Distance Box plot

{#if mounted}
<BoxHullDistErrors />
{/if}

> @label:fig:box-hull-dist-errors Box plot of interquartile ranges (IQR) of hull distance errors for each model. The whiskers extend to the 5th and 95th percentiles. The horizontal line inside the box shows the median. BOWSR has the highest median error, while Voronoi RF has the highest IQR. Note that MEGNet and CGCNN are the only models with a positive median. Their hull distance errors are biased towards more frequently predicting thermodynamic instability, explaining why they are closest to getting the overall number of stable structures in the test set right (see cumulative precision/recall in @fig:cumulative-precision-recall).

BOWSR has the largest median error, while Voronoi RF has the largest IQR. Note that MEGNet and CGCNN are the only models with a positive median. Their hull distance errors are biased towards more frequently predicting thermodynamic instability, explaining why they are closest to correctly predicting the overall number of stable structures in the test set (see @fig:cumulative-precision-recall).

### Classification Histograms using Model-Predicted Energies

{#if mounted}
<HistClfPredHullDistModels />
{/if}

> @label:fig:hist-clf-pred-hull-dist-models Distribution of model-predicted hull distance colored by stability classification. Models are sorted from top to bottom by F1 score. The thickness of the red and yellow bands shows how often models misclassify as a function of how far away from the convex hull they place a material. While CHGNet and M3GNet perform almost equally well overall, these plots reveal that they do so via different trade-offs. M3GNet commits fewer false negatives but more false positives predictions compared to CHGNet. In a real discovery campaign, false positives have a higher opportunity cost than false negatives, since they result in wasted DFT relaxations or even synthesis time in the lab. A false negative by contrast is just one missed opportunity out of many. This observation is also reflected in the higher TPR and lower TNR of M3GNet vs CHGNet in @tab:metrics-table-uniq-protos, as well as the lower rolling MAE for CHGNet vs M3GNet on the stable side (left half) of @fig:rolling-mae-vs-hull-dist-models and vice-versa on the unstable side (right half).

Note the CGCNN+P histogram is more strongly peaked than CGCNN's which agrees better with the actual DFT ground truth [distribution of hull distances](/data#--target-distribution) in our test set. This explains why CGCNN+P performs better as a regressor, but also reveals how it can perform simultaneously worse as a classifier. By moving predictions closer to the stability threshold at 0 eV/atom above the hull, even small errors are significant enough to tip a classification over the threshold.

### Measuring extrapolation performance from WBM batch robustness

As a reminder, the WBM test set was generated in 5 successive batches, each step applying another element replacement to an MP source structure or a new stable crystal generated in one of the previous replacement rounds. The likelihood by which one element replaces another is governed by ISCD-mined chemical similarity scores for each pair of elements. Naively, one would expect model performance to degrade with increasing batch count, as repeated substitutions should on average 'diffuse' deeper into uncharted regions of material space, requiring the model to extrapolate more. We observe this effect for some models much more than others.

@Fig:rolling-mae-vs-hull-dist-wbm-batches-models shows the rolling MAE as a function of distance to the convex hull for each of the 5 WBM rounds of elemental substitution. These plots show a stronger performance decrease for Wrenformer and Voronoi RF than for UIPs like M3GNet, CHGNet, MACE and even force-less GNNs with larger errors like MEGNet and CGCNN.

{#if mounted}
{@const style=`aspect-ratio: 1.3;`}

  <div style="display: grid; grid-template-columns: 1fr 1fr; margin: 0 -1em 0 -4em;">
    <RollingMaeVsHullDistWbmBatchesOrb {style} />
    <RollingMaeVsHullDistWbmBatchesMatterSim {style} />
    <RollingMaeVsHullDistWbmBatchesGnome {style} />
    <RollingMaeVsHullDistWbmBatchesSevenNet {style} />
    <RollingMaeVsHullDistWbmBatchesMace {style} />
    <RollingMaeVsHullDistWbmBatchesChgnet {style} />
    <RollingMaeVsHullDistWbmBatchesM3gnet {style} />
    <RollingMaeVsHullDistWbmBatchesAlignn {style} />
    <RollingMaeVsHullDistWbmBatchesMegnet {style} />
    <RollingMaeVsHullDistWbmBatchesWrenformer {style} />
    <RollingMaeVsHullDistWbmBatchesCgcnn {style} />
    <RollingMaeVsHullDistWbmBatchesVoronoiRf {style} />
  </div>
{/if}

> @label:fig:rolling-mae-vs-hull-dist-wbm-batches-models Rolling MAE as a function of distance to the convex hull for different models.
> On WBM batch 1, Wrenformer performs comparably to the SOTA UIPs M3GNet and CHGNet.
> However, the M3GNet and CHGNet UIPs show stronger extrapolative performance, as they show much smaller deterioration in performance on later batches that move further away from the original MP training distribution.
> Wrenformer, by contrast, exhibits a pronounced increase in MAE with batch count.
> We view these plots as a strong indicator that Matbench Discovery is indeed testing out-of-distribution extrapolation performance as it is Occam's
> razor explanation for the observed model performance drop with increasing batch count.
> The effect is larger for some models than others but batches 4 and 5 consistently incur the highest convex hull distance errors for all models.

MEGNet and CGCNN both incur a higher rolling MAE than Wrenformer across all 5 batches and across most or all of the hull distance range visible in these plots. However, similar to the UIPs, MEGNet and CGCNN exhibit very little degradation for higher batch counts. The fact that higher errors for later batches are specific to Wrenformer and Voronoi RF suggests that training only on composition and coarse-grained structural features (spacegroup and Wyckoff positions in the case of Wrenformer; coordination numbers, local environment properties, etc. in the case of Voronoi RF) alone is insufficient to learn an extrapolatable map of the PES.

Given its strong performance on batch 1, it is possible that given sufficiently diverse training data, Wrenformer could become similarly accurate to the UIPs across the whole PES landscape at substantially less training and inference cost. However, the loss of predictive forces and stresses may make Wrenformer unattractive for certain applications even then.

### Largest Errors vs. DFT Hull Distance

Given the large variety of models tested, we asked whether any additional insight into the errors can be gained from looking at how the predictions vary between different models.
In @fig:scatter-largest-errors-models-mean-vs-true-hull-dist we see two distinct groupings emerge when looking at the 200 structures with the largest errors.
This clustering is particularly apparent when points are colored by model disagreement.

{#if mounted}
<ScatterLargestErrorsModelsMeanVsTrueHullDist />
{/if}

> @label:fig:scatter-largest-errors-models-mean-vs-true-hull-dist DFT vs predicted hull distance (average over all models) for the 200 largest error structures.
> Points are colored by model disagreement as measured by the standard deviation in hull distance predictions from different models.
> Point scale with the number of atoms in the structures.
> This plot shows that high-error predictions are biased towards predicting too small hull distances.
> This is unsurprising considering MP training data mainly consists of low-energy structures.
> There is a strong color separation between the mostly dark blue low-energy bias predictions and light blue/green high-error predictions.
> Blue means models are in good agreement, i.e. all models are "wrong" together.
> Green means large-error predictions with little model agreement, i.e. all models are wrong in different ways.
> Some of the blue points with large errors yet good agreement among models may be accurate ML predictions for a DFT relaxation gone wrong.
> The dark blue points also tend to be larger corresponding to larger structures where DFT failures are less surprising.
> This suggests ML model committees might be used to cheaply screen large databases for DFT errors in a high-throughput manner.

### MEGNet formation energies from UIP-relaxed structures

{#if mounted}
<MetricsTableMegnetUipCombos />
{/if}

> @label:tab:metrics-table-megnet-uip-combos Except for MEGNet RS2RE, all rows in this table show metrics for the task of predicting the relaxed energy given the unrelaxed structure (IS2RE). Metrics in rows labeled M3GNet→MEGNet and CHGNet→MEGNet are the result of passing M3GNet/CHGNet-relaxed structures into MEGNet for formation energy prediction. Both model combos perform worse than using the respective UIPs on their own with a more pronounced performance drop from CHGNet to CHGNet→MEGNet than M3GNet to M3GNet→MEGnet. This suggests that MEGNet has no additional knowledge of the PES that is not already encoded in the UIPs. However, both combos perform better than MEGNet on its own, demonstrating that ML relaxation imparts some of the latent knowledge of the PES into the relaxed structure from which a more informed final energy prediction is possible. This underscores the utility of ML relaxation at a very low cost for any downstream structure-dependent analysis.

Interestingly, @tab:metrics-table-megnet-uip-combos shows that MEGNet RS2RE performs worse than CHGNet and M3GNet even when it has access to the ground truth DFT relaxed structures when making relaxed energy predictions.
While the regression and classification metrics for MEGNet RS2RE are better than M3GNet→MEGNet and CHGNet→MEGNet, the difference is surprisingly small, especially compared to the regression performance of standalone MEGNet.
This suggests that for downstream use, at least by an imprecise method such as another ML model, the UIP-relaxed structures are almost as helpful as the DFT-relaxed ones.

CHGNet, M3GNet and MEGNet were all trained on slightly different targets.
While MEGNet was trained to predict formation energies, M3GNet and CHGNet were trained on raw DFT energies with the important distinction that CHGNet targets include [MP2020 energy correction scheme](https://github.com/materialsproject/pymatgen/blob/02a4ca8aa0277b5f6db11f4de4fdbba129de70a5/pymatgen/entries/compatibility.py#L823) [@wang_framework_2021] while M3GNet targets do not.

Changing the target from raw DFT energies to formation energies in principle merely constitutes a change of gauge since the difference is just a linear transformation of subtracting the elemental reference energies weighted by composition.
But in practice, one might expect the problem of stability prediction to be easier when trained on formation energies, as the model output is one less step removed from the final target of convex hull distance.
However, empirical evidence so far suggests the opposite.
Both UIP→MEGNet combos perform worse than those same UIPs by themselves.
There are too many confounding effects at play to draw firm conclusions but this tentatively suggests using formation energy as targets offers no benefits over raw DFT.

### Exploratory Data Analysis

<img src={wbm_matminer_2d_umap_url} alt="WBM matminer features 2D UMAP projection" style="background-color: rgba(255, 255, 255, 0.3); border-radius: 4px; max-width: 550px; margin: auto; display: block;">

> @label:fig:wbm-final-struct-matminer-features-2d-umap 2D UMAP projection of `matminer` [@ward_matminer_2018] features computed for all initial structures in the WBM test set.
> Points are colored by elemental substitution step.
> Step 1 means a single element replacement applied to a source structure taken from MP.
> Step 2 means two iterations, i.e. two successive element replacements and so on.
> Step 0 are relaxed MP structures.
> While the MP structures are not disjoint from WBM, there are large areas of the WBM UMAP projection that are not represented in MP, indicating our benchmark requires some degree of out-of-distribution extrapolation to perform well.

To give high-level insights into the MP training and WBM test set used in this work, we include element distributions for structures in both data sets (@fig:element-counts-by-occurrence).
To show how frame selection from MP structure relaxation affected relative elemental abundance between MP relaxed structures and the snapshots in MPtrj, @fig:element-counts-ratio-by-occurrence shows element occurrence ratios between MPtrj and MP.
Believing MPtrj to be an influential dataset for the near-term continued development of universal interatomic potentials, we plot histograms showing the distributions of target values for energies, forces, stresses and magnetic moments in @fig:mp-trj-hists.
Similarly, @fig:wbm-energy-hists plots the distribution of formation energies and convex hull distances (with respect to the convex hull spanned by MP materials only) of the WBM test set.
This plot not only gives insight into the nature of the dataset but also emphasizes the increased difficulty of stability vs formation energy prediction arising from the much narrower distribution of convex hull distances compared to the more spread-out formation energy distribution.

To get a visual sense of coverage of atomic configuration space afforded by WBM compared to MP (as an indicator for the degree of out-of-distribution extrapolation required), we project `matminer` [@ward_matminer_2018] features for all initial structures as well as relaxed MP structures onto two UMAP dimensions in @fig:wbm-final-struct-matminer-features-2d-umap.

Finally, @fig:mp-vs-mp-trj-vs-wbm-arity-hist shows the elements-per-structure distribution of MP, MPtrj and WBM, normalized by dataset size.
The mode of all three datasets is 3, but WBM's share of ternary phases is noticeably more peaked than MP's, which includes small numbers of unary and senary phases.
The histogram in @fig:mp-trj-n-sites-hist shows the distribution of the number of sites in MPtrj structures. The inset displays the same histogram log-scaled y-axis as well as a cumulative line to show that 90% of MPtrj structures contain fewer than 70 sites.

{#if mounted}

<!-- <ElementCountsByOccurrence /> -->

{/if}

> @label:fig:element-counts-by-occurrence The number of structures containing a given element in the MP training set, the MPtrj dataset containing multiple from every relaxation trajectory in MP [@deng_chgnet_2023], and WBM test set [@wang_predicting_2021].
> The WBM test set is noticeably more chemically diverse than MP (and, by extension, MPtrj).
> Made with pymatviz [@riebesell_pymatviz_2022].

{#if mounted}

<!-- <ElementCountsRatioByOccurrence /> -->

{/if}

> @label:fig:element-counts-ratio-by-occurrence **)** shows the ratio of elements in the WBM test set to the MP training set.
> Similarly, **)** shows the ratio of elements in the MPtrj dataset to the MP training set.
> We note a slight overabundance of structures containing hydrogen and halides, indicating that more frames were selected from structures containing these elements which might correlate with the number of ionic steps to find their ground states.

{#if mounted}
<MPvsMPTrjVsWbmArityHist />
{/if}

> @label:fig:mp-vs-mp-trj-vs-wbm-arity-hist Distribution of unique elements per structure in MP, MPtrj and WBM.
> The bar heights are normalized by the total number of structures in each data set.
> WBM is dominated by ternary phases making up 74% of the data set followed by about 13% each of binaries and quaternaries.
> MP has a more even distribution, in particular with more than double the relative share of quaternary phases and a significant number of quinternaries which are almost absent from WBM.
> Not shown in this plot for visual clarity are 3% of MP structures containing more than 5 elements (up to 9).
> We also include MPtrj in this plot to show a slight drop in the relative abundance of quinternaries and higher phases vs MP ground states.
> This may be due to a poor choice of convergence criteria in early MP relaxation workflows that scaled with the size of the structure (see `EDIFF_PER_ATOM` parameter in `pymatgen` VASP input sets), resulting in unconverged large structures with short relaxation trajectories entering the database.
> Short relaxations would result in fewer frames of such structures selected for MPtrj.
> This assumes structures of higher arity correlate with larger structures.

{#if mounted}
<MpTrjNSitesHist />
{/if}

> @label:fig:mp-trj-n-sites-hist Histogram of number of atoms per structure.
> The inset shows the same distribution log-scaled to visualize the tail of large structures.
> The green cumulative line in the inset shows that 82% have less than 50 sites and 97% of structures in MPtrj have less than 100 atoms.

{#if mounted}
{@const style=`aspect-ratio: 3/2;`}

  <div style="display: grid; gap: 1em; grid-template-columns: 1fr 1fr;">
    <MPTrjEFormHist {style} />
    <MPTrjForcesHist {style} />
    <MPTrjStressesHist {style} />
    <MPTrjMagMomsHist {style} />
  </div>
{/if}

> @label:fig:mp-trj-hists Distribution of energies, forces, stresses and magnetic moments MPtrj.
> The bimodality in the formation energy distribution is due to the MP anion correction scheme [@wang_framework_2021] which significantly lowers some formation energies, especially for oxides.

{#if mounted}
{@const style=`aspect-ratio: 1.2;`}

  <div style="display: grid; gap: 1em; grid-template-columns: 1fr 1fr;">
    <HistWbmEForm {style} />
    <HistWbmHullDist {style} />
  </div>
{/if}

> @label:fig:wbm-energy-hists Distribution of formation energies and convex hull distances (both per atom) in the WBM test set [@wang_predicting_2021].
> **Left**: We removed what we believe to be invalid calculations with formation energies outside the range of ±5 eV/atom plotted in.
> This filter discarded 22 structures to the right of the visible distribution with formation energies between 5 eV/atom and 80 eV/atom as well as 502 suspicious structures to the left of it with exactly −10 eV/atom formation energy.
> However, even the cleaned formation energy distribution still exhibits much wider spread than the convex hull distance distribution (**right**), spanning almost 10 eV/atom vs less than 1 eV/atom spread in the hull distances.
> This highlights why stability prediction is a much more challenging task than predicting energy of formation.
> It requires correctly ranking the subtle energy differences between chemically similar compounds in the same chemical system rather than comparing a single material with the reference energies of its constituent elements.
> DFT has been shown to significantly benefit from the systematic cancellation of errors between chemically similar systems when trying to identify the lowest-lying polymorph [@hautier_accuracy_2012].
> This beneficial cancellation has yet to be conclusively demonstrated for ML stability predictions.
> So far, only the lack thereof has been shown in [@bartel_critical_2020] where they encountered a much more random error distribution among similar chemistries than simulations from first principles.

<!-- {#if mounted}
<MpTrjPtableHists />
{/if}

> @label:fig:mp-trj-ptable-hists Distribution of magnetic moments and forces for each element MPtrj. This data is used as training targets for all interatomic potentials in this work (only CHGNet uses the absolute value of magnetic moments as targets).
> The number in the top right corner of each element tile counts the number of target values for that element in all of MPtrj.
> $y$-axes are log-scaled to reveal the tail of high magnetic moments in some elements.
> ) reveals rare erroneous data points in MPtrj.
> For instance, has a single-point calculation with a highly unphysical magnetic moment of 17$\mu_\text{B}$.
> For visualization purposes, the $y$-axes are again log-scaled and distributions are truncated at 10 eV/Å.
> Oxygen has the largest outliers with mean absolute forces of up to 160 eV/Å. -->

### Predicted Formation Energy Analysis

We show the distribution for all materials in the WBM test set of model-predicted formation energies vs DFT ground truth in @fig:e-form-parity-models

<!-- and the distribution of convex hull distance errors projected onto elements by composition for MACE and CHGNet in @fig:ptable-each-error-hists. -->

{#if mounted}
<EFormParityModels />
{/if}

> @label:fig:e-form-parity-models Parity plots of model-predicted formation energies vs DFT ground truth, color-coded by log density of points.
> The models are sorted from left to right and top to bottom by MAE.
> While similar to the parity plots in @fig:each-parity-models which shows the predicted convex hull distance vs PBE hull distance, this figure better visualizes the point density due to formation energy's wider spread.
> We observe broadly the same failure modes with occasional high DFT energy outliers predicted as near 0 formation energy by the models.

<!-- {#if mounted}
<PtableEachErrorHists />
{/if}

> @label:fig:ptable-each-error-hists Distribution of convex hull distance errors for all materials in the WBM test set projected onto elements by composition, for the two best models we tested, MACE [@batatia_mace_2023] and CHGNet [@deng_chgnet_2023].
> Element error projection is best explained by example: assume has a convex hull distance error of 50 meV/atom, $\frac{2}{5}$ or 20 meV/atom of that error would be attributed to and $\frac{3}{5}$ or 30 meV/atom to .
> Thus would show up as one count of 20 for and one count of 30 for in this histogram.
> We observe very similar error distributions for MACE and CHGNet, confirming that the two models are learning a similar map of the PES.
> There's a weak trend of distributions more sharply peaked at 0 error for CHGNet than MACE.
> We also observe slightly skewed error distributions from MACE on the halides which accumulate more density on negative than positive errors.
> That is, MACE tends to underestimate the convex hull distance and therefore overestimate the stability of halide-containing materials. -->
