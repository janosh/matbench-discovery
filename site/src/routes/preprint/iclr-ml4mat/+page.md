<script>
  import { repository as repo } from '$site/package.json'
  import MetricsTable from '$figs/metrics-table.svelte'
  import CumulativePrecisionRecall from '$figs/cumulative-precision-recall.svelte'
  import RollingMaeVsHullDistModels from '$figs/rolling-mae-vs-hull-dist-models.svelte'
  import { browser } from '$app/environment'
</script>

> This is a 4-page extended abstract submitted to [ICLR ML4Materials workshop](https://ml4materials.com) on 2023-02-10.

<summary>

We present a new machine learning (ML) benchmark for materials stability predictions named `Matbench Discovery`. A goal of this benchmark is to highlight the need to focus on metrics that directly measure their utility in prospective discovery campaigns as opposed to analyzing models based on predictive accuracy alone. Our benchmark consists of a task designed to closely simulate the deployment of ML energy models in a high-throughput search for stable inorganic crystals. We explore a wide variety of models covering multiple methodologies ranging from random forests to GNNs, and from one-shot predictors to iterative Bayesian optimizers and interatomic potential relaxers. We find M3GNet to achieve the highest F1 score of 0.58 and $R^2$ of 0.59 while MEGNet wins on discovery acceleration factor (DAF) with 2.94. Our results provide valuable insights for maintainers of high throughput materials databases to start using these models as triaging steps to more effectively allocate compute for DFT relaxations.

</summary>

## Introduction

For nearly two decades, ever since the work of Behler and Parrinello @behler_generalized_2007 which introduced a custom neural network for learning the density-functional theory (DFT) potential energy surface (PES), material scientists have devoted significant effort to developing custom model architectures for tackling the problem of learning the PES.
Initially, most of these models were trained and deployed as interatomic potentials to study known materials of interest which required curating custom training data for each application @bartok_machine_2018 @deringer_general-purpose_2020.
As larger and more diverse datasets emerged from initiatives like the Materials Project (MP) @jain_commentary_2013 or the Open Quantum Materials Database (OQMD) @saal_materials_2013, researchers have begun to train models that cover the full periodic table opening up the prospect of ML-guided materials discovery.

Yet despite many advances in ML for materials discovery, it is unclear which methodology performs best at predicting material stability.
Recent areas of progress include

1. one-shot predictors like Wren @goodall_rapid_2022,
1. universal force predictors such as M3GNet @chen_universal_2022 that emulate density functional theory to relax crystal structures according to Newton's laws, and
1. Bayesian optimizers like BOWSR that, paired with any energy model, treat structure relaxation as a black-box optimization problem @zuo_accelerating_2021.

In this work, we aim to answer which of these is the winning methodology in a future-proof benchmark that closely simulates using ML to guide a real-world discovery campaign.

## Related Work

### Using Stability rather than Formation Energies

In 2020, Chris Bartel et al. @bartel_critical_2020 benchmarked 7 models, finding all of them able to predict DFT formation energies with useful accuracy.
However, when asked to predict stability (specifically decomposition enthalpy), the performance of all models deteriorated sharply.
This insight meant that ML models are much less useful than DFT for discovering new solids than prior studies had suggested.
The paper identified two main reasons for the sharp decline in predictive power:

1. Stability is a property not only of the material itself but also the chemical space of competing phases it occupies. Current ML algorithms are given an input that only describes the single material they are asked to predict, leaving them clueless of competing phases.
1. Unlike DFT, ML models appear to benefit less from systematic error cancellation across similar chemistries.

Bartel et al. showed that to demonstrate the utility of ML for materials discovery, the vanity metric of formation energy accuracy must be replaced with stability predictions.
Moreover, the qualitative leap in performance from Roost @goodall_predicting_2020, the best compositional model benchmarked, to CGCNN @xie_crystal_2018, the single structural model they tested, shows structure plays a crucial role in determining the stability of materials.
However, using the DFT-relaxed structure as input to CGCNN renders the discovery pipeline circular as the input becomes the very thing we aim to find. A true test of prospective utility requires using unrelaxed structures as the next most information-rich input.

### Matbench

As the name suggests, this work seeks to expand upon the original Matbench suite of property prediction tasks @dunn_benchmarking_2020. By providing a standardized collection of datasets along with canonical cross-validation splits for model evaluation, Matbench helped focus the field of ML for materials, increase comparability across papers
and attempt to accelerate the field similar to what ImageNet did for computer vision.

Matbench released a test suite of 13 supervised tasks for different material properties ranging from thermal (formation energy, phonon frequency peak), electronic (band gap), optical (refractive index) to tensile and elastic (bulk and shear moduli).
They range in size from ~300 to ~132,000 samples and include both DFT and experimental data sources.
Importantly, all tasks were exclusively concerned with the properties of known materials.
We believe a task that simulates a materials discovery campaign by requiring materials stability predictions from unrelaxed structures to be a missing piece here.

### The Open Catalyst Project

The Open Catalyst Project (OCP) is a large-scale initiative to discover substrate-absorbate combinations that catalyze key industrial reactions processing said absorbates into more useful products.
The OCP has released two data sets thus far, OCP20 @chanussot_open_2021 and OCP22 @tran_open_2022, that can be used for training and benchmarking ML models.

However, the ambition and scale of OCP comes with limitations for its use as a benchmark.
OCP20 is already 10x larger than the largest crystal structure data sets available imposing a high barrier to entry for researchers without access to cloud-scale computing.
In contrast, we believe the discovery of stable materials is a problem where ML methods have matured enough to be usefully deployed at scale after training for only $\mathcal{O}(10^2)$ GPU hours.

## Data Sets

The choice of data for the train and test sets of this benchmark fell on the latest Materials Project (MP) @jain_commentary_2013 database release ([2022.10.28] at time of writing) and the WBM dataset @wang_predicting_2021.

[2022.10.28]: https://docs.materialsproject.org/changes/database-versions#v2022.10.28

### The Materials Project - Training Set

The Materials Project is a well-known effort to calculate the properties of all inorganic materials using high-throughput ab-initio methods.
At the time of access, the Materials Project database contains approximately 154k crystals (providing relaxed+initial structure and the relaxation trajectory for each of them) covering a diverse range of chemistries.
For our benchmark, the training set is all data available from the [2022.10.28] MP release. Models are free to train on relaxed and/or unrelaxed structures or the full DFT relaxation trajectory. This flexibility is intended to allow authors to experiment and exploit the large variety of data available.

### WBM - Test Set

The WBM data set @wang_predicting_2021 consists of ~257k structures generated via chemical similarity-based elemental substitution of MP source structures followed by DFT relaxation and convex hull distance calculation.
Throughout this work, we define stability in terms of being on or below the convex hull of the MP training set. ~42k out of ~257k materials in WBM satisfy this criterion.
As WBM explores regions of materials space not well sampled by MP, many of these materials discovered that are stable w.r.t. MP's convex hull are not stable with respect to each other.
Only around ~20k were found to remain on the convex hull when merging the MP and WBM hulls.
This observation highlights a critical aspect of this benchmark in that we purposely operate with an incomplete convex hull. Only current knowledge is accessible to a real discovery campaign. Hence our metrics are designed to reflect this.

Moreover, to simulate a discovery campaign our test set inputs are unrelaxed structures obtained from chemical substitution of MP source structures but our target labels are the relaxed PBE formation energies. This opens up the opportunity to explore how different approaches (one-shot, force-based pseudo-relaxation, black-box pseudo-relaxation, etc.) compare for materials discovery.

## Models

Our initial benchmark release includes 8 models.

1. **Voronoi+RF** @ward_including_2017 - A random forest trained to map a combination of composition-based Magpie features and structure-based relaxation-invariant Voronoi tessellation features (effective coordination numbers, structural heterogeneity, local environment properties, ...) to DFT formation energies.

1. **Wrenformer** @goodall_rapid_2022 - For this benchmark, we introduce Wrenformer which is a variation on Wren @goodall_rapid_2022 constructed using standard QKV-Transformer blocks to reduce memory usage, allowing it to scale to structures with >16 Wyckoff positions.

1. **CGCNN** @xie_crystal_2018 - The Crystal Graph Convolutional Neural Network (CGCNN) was the first neural network model to directly learn 8 different DFT-computed material properties from a graph representing the atoms and bonds in a periodic crystal.

1. **CGCNN+P** @gibson_data-augmentation_2022 - This work proposes a simple, physically motivated structure perturbations to augment stock CGCNN's training data of relaxed structures with structures resembling unrelaxed ones but mapped to the same DFT final energy. Here we chose $P=5$, meaning the training set was augmented with 5 random perturbations of each relaxed MP structure mapped to the same target energy.

1. **MEGNet** @chen_graph_2019 - MatErials Graph Network is another GNN similar to CGCNN for material properties of relaxed structures that also updates the edge and global features (like pressure, temperature, entropy) in its message passing operation.

1. **M3GNet** @chen_universal_2022 - M3GNet is a GNN-based universal (as in full periodic table) interatomic potential (UIP) for materials trained on up to 3-body interactions in the initial, middle and final frame of MP DFT relaxations. The model takes the unrelaxed input and emulates structure relaxation before predicting energy for the pseudo-relaxed structure.

1. **BOSWR + MEGNet** @zuo_accelerating_2021 - BOWSR uses a symmetry-constrained Bayesian optimizer (BO) with a surrogate energy model (here MEGNet) to perform an iterative exploration-exploitation-based search of the potential energy landscape. The high sample count needed to explore the PES with BO makes this by far the most expensive model.

## Results

<MetricsTable />

> @label:fig:metrics-table Regression and classification metrics for all models tested on our benchmark. The heat map ranges from yellow (best) to blue (worst) performance. DAF = discovery acceleration factor (see text), TPR = true positive rate, TNR = true negative rate, MAE = mean absolute error, RMSE = root mean squared error

@Fig:metrics-table shows performance metrics for all models included in the initial release of Matbench Discovery.
M3GNet takes the top spot on most metrics and emerges as current SOTA for ML-guided materials discovery. The discovery acceleration factor (DAF) measures how many more stable structures a model found among the ones it predicted stable compared to the dummy discovery rate of 43k / 257k $\approx$ 16.7% achieved by randomly selecting test set crystals. Consequently, the maximum possible DAF is ~6. This highlights the fact that our benchmark is made more challenging by deploying models on an already enriched space with a much higher fraction of stable structures over randomly exploring materials space. As the convex hull becomes more thoroughly sampled by future discovery, the fraction of unknown stable structures decreases, naturally leading to less enriched future test sets which will allow for higher maximum DAFs. The reason MEGNet outperforms M3GNet on DAF becomes clear from @fig:cumulative-precision-recall by noting that MEGNet's line ends closest to the total number of stable materials. The other models overpredict this number, resulting in large numbers of false positive predictions that drag down their DAFs.

{#if browser}
<RollingMaeVsHullDistModels />
{/if}

> @label:fig:rolling-mae-vs-hull-dist-models Rolling MAE on the WBM test set as the energy to the convex hull of the MP training set is varied. The white box in the bottom left indicates the size of the rolling window. The highlighted 'triangle of peril' shows where the models are most likely to misclassify structures. We only show the 6 best performing models for visual clarity. We only show the 6 best performing models for visual clarity.

{#if browser}
<CumulativePrecisionRecall style="margin: 0 -2em 0 -4em;" />
{/if}

> @label:fig:cumulative-precision-recall Cumulative precision and recall over the course of a simulated discovery campaign. This figure highlights how different models will perform best depending on the setup of the screening campaign.

@Fig:rolling-mae-vs-hull-dist-models visualizes a model's reliability as a function of a material's hull distance. The lower its rolling MAE exits the shaded triangle, the better. Inside this area, the model's mean error is larger than the distance to the convex hull, making misclassifications likely. Outside the triangle even if the model's error points toward the stability threshold at 0 eV from the hull (the plot's center), the mean error is too small to move a material over the stability threshold which would cause a false stability classification. M3GNet achieves the lowest overall MAE and exits the peril zone much sooner than other models on the right half of the plot. This means it rarely misclassifies unstable materials that lie more than 40 meV above the hull. On the plot's left half, CGCNN+P exits the peril zone first, albeit much further from the hull at more than 100 meV below. Essentially, all models are prone to false negative predictions even for materials far below the known hull. We note that while F1 score and DAF of models that make one-shot predictions directly from unrelaxed inputs (CGCNN & MEGNet) are seemingly unaffected, the $R^2$ of these models is significantly worse.

Despite their low accuracy in one-shot predicting relaxed energies from unrelaxed structures, CGCNN and MEGNet achieve high F1 scores and DAFs. This can be explained by unrelaxed inputs being in higher energy configurations than their relaxed counterparts which makes them more likely to be unstable with respect to the training set convex hull, even in cases where the relaxed structure is stable. This biases one-shot GNNs towards predicting unrelaxed input structures as unstable, resulting in higher true negative rates that offset the lower true positive rates in the F1 and DAF metrics. This phenomenon is expected due to the training and testing configuration mismatch and explains the success of previous screening attempts that used such GNN models for screening unrelaxed inputs @park_developing_2020.

## Discussion

From @fig:metrics-table we see several models achieve a DAF > 2 in this realistic benchmark scenario.
Consequently, the benefits of deploying ML-based triage in high-throughput computational materials discovery applications likely warrant the time and setup required.
However, there are many aspects on which further progress is necessary, for example, models still make large numbers of false positive predictions for materials over 50 meV above the convex hull and much less likely to be synthesizable, greatly reducing the DAF.
The results obtained from version 1 of our benchmark show that ML universal interatomic potentials like M3GNet are the most promising methodology to pursue going forward, being both ~20x cheaper to run than black box optimizers like BOWSR and having access to more training structures than coordinate-free approaches like Wrenformer.

Although the task of discovery will necessarily become more challenging over time as currently undersampled regions of materials space are explored, the path to making ML a ubiquitous discovery tool appears straightforward and is one the field is already pursuing: training foundational UIPs on significantly more data may get us there even without further algorithmic or model improvements.
We welcome further model submissions as well as data contributions for version 2 of this benchmark to the GitHub repo at
[{repo}]({repo}).

## Acknowledgments

Janosh Riebesell acknowledges support from the German Academic Scholarship Foundation ([Studienstiftung](https://wikipedia.org/wiki/Studienstiftung)).

A big thank you to

- Hai-Chen Wang and co-authors for creating and freely providing the WBM data set
- Jason Blake Gibson, Shyue Ping Ong, Chi Chen, Tian Xie, Bowen Deng, Peichen Zhong, Ekin Dogus Cubuk for helpful discussions
- Philipp Benner ([@pbenner](https://github.com/pbenner)) for [finding and reporting many bugs]({repo}/issues?q=is%3Aissue+author%3Apbenner+) in the data loading routines before the v1 release.

## Author Contributions

Janosh Riebesell: Methodology, Software, Data Curation, Training and Testing Models, Formal Analysis. Rhys Goodall: Conceptualization, Software, Formal Analysis. Anubhav Jain: Supervision. Kristin Persson: Supervision. Alpha Lee: Supervision.
