<script>
  import { repository as repo } from '$site/package.json'
  import EachScatter from '$figs/each-scatter-models.svelte'
  import MetricsTable from '$figs/metrics-table.svelte'
  import CumulativeClfMetrics from '$figs/cumulative-clf-metrics.svelte'
  import RollingMaeModels from '$figs/rolling-mae-vs-hull-dist-models.svelte'
  import ElementErrorsPtableHeatmap from '$models/element-errors-ptable-heatmap.svelte'
  import { browser } from '$app/environment'
</script>

<summary>

We present a new machine learning (ML) benchmark for materials stability predictions named `Matbench Discovery`. A goal of this benchmark is to highlight the need to focus on metrics that directly measure their hit rate in prospective discovery campaigns as opposed to analyzing models based on predictive accuracy alone. Our benchmark consists of a task designed to closely simulate the deployment of ML energy models in a high-throughput search for stable inorganic crystals. To shed light on the question which type of ML performs best at materials discovery, we explore a wide variety of models covering multiple methodologies. Our selection ranges from random forests to GNNs, from one-shot predictors over iterative Bayesian optimizers to interatomic potential (IAP) relaxers that closely emulate DFT. We find the M3GNet IAP to achieve the highest F1 score of 0.58 and $R^2$ of 0.59 while MEGNet wins on discovery acceleration factor (DAF) with 2.94. Our results provide valuable insights for maintainers of high throughput materials databases to start using these models as triaging steps for effectively allocating DFT relaxations.

</summary>

## Introduction

### Why Materials Discovery?

@davies_computational_2016 identified $~10^{10}$ possible quaternary materials using electronegativity and charge-based allowed by electronegativity and charge-balancing rules. Of these, only $~10^5$ are known experimentally and $~10^6$ have been simulated. The space of quinternaries and higher is even less explored, leaving vast numbers of potentially useful materials to be discovered. The discovery of new materials is a key driver of technological advancement and lies on the critical path to more efficient solar cells, lighter and longer-lived batteries, and smaller more efficient transistor gates needed for further device miniaturization. The materials we have mastered define our technological abilities so much so that we name eras of human development after them, going from stone and bronze age via the iron age to our current silicon age. In light of global warming, these advances cannot come fast enough.

### Why use machine learning for discovery?

Despite significant advances in empirical, theoretical and computational materials science, discovering new materials still requires labor-intensive experimentation, complex calculations, trial-and-error and often happens fortuitously rather than through rational design. ML methods can be useful for materials discovery because they efficiently extract and distill trends from very large datasets, can handle high dimensionality, multiple objectives, uncertainty, and noisy or sparse data.

Statistical models are also orders of magnitude faster than ab-initio simulation. As such, they are most suitable for high-throughput searches as a pre-filter to more expensive, higher-fidelity simulation methods. The use of neural networks for learning the density-functional theory (DFT) potential energy surface (PES) can be traced as far back as @behler_generalized_2007. This work opened the floodgates for material scientists to devote significant effort into fitting ever more sophisticated models to the PES.
Initially, most of these models were trained and deployed as interatomic potentials to study known materials of interest which required curating custom training data for each application @bartok_machine_2018 @deringer_general-purpose_2020.
As larger and more diverse datasets emerged from initiatives like the Materials Project (MP) @jain_commentary_2013 or the Open Quantum Materials Database (OQMD) @saal_materials_2013, researchers have begun to train models that cover the full periodic table, opening up the prospect of ML-guided materials discovery to increase hit rate and speed of DFT/expert-driven searches.

### Limitations of current benchmarking in ML

Yet despite many advances in ML for materials it is unclear which methodology performs best at predicting material stability, let alone which model. This is due to the lack of a standardized benchmark task that accurately simulates applying models in a prospective materials discovery campaign.

1. **Lack of realism**: Benchmark tasks can be idealized and simulate overly simplified conditions that do not reflect the real-world challenges a model is expected to overcome when used in an actual discovery campaign. This can lead to pretty leaderboards listing seemingly SOTA models that underwhelm when used in production. Examples of how this comes about are choosing the wrong target or picking an unrepresentative train/test split.
1. **Limited diversity**: Benchmark datasets may be too small and contain only a limited number of materials, unrepresentative of the huge diversity of materials space. This can make models look good even if they fail to generalize.
1. **Opportunity cost**: Bad benchmarks give insufficient consideration to the cost of a failed experiment. Looking purely at global metrics like $\text{MAE}$, $\text{RMSE}$ and $R^2$ can give practitioners a false sense of security. Precision, recall F1 It has been shown that even accurate models are susceptible to unexpectedly high false-positive rates that can cause experimentalists to waste their time and resources. Many benchmark tasks do not consider the cost or practicality of synthesizing the materials, which is an important aspect in the discovery of new materials.
1. **Scalability**: Many benchmark tasks have too small data to adequately simulate the high-throughput and large-date regimes that future discovery efforts are likely to encounter. Confining model testing to the small data regime can obfuscate poor scaling relations like Gaussian Processes whose training costs grow cubically with training sample count or random forests that achieve outstanding performance on few data points but fail to extract the full information content out of larger datasets, leading to flatter learning curves compared to neural networks and worse performance when large amounts of training data are available.

Recent areas of progress include

1. one-shot predictors like Wren @goodall_rapid_2022,
1. universal force predictors such as M3GNet @chen_universal_2022 that emulate density functional theory to relax crystal structures according to Newton's laws, and
1. Bayesian optimizers like BOWSR that, paired with any energy model, treat structure relaxation as a black-box optimization problem @zuo_accelerating_2021.

Ideally, the question of which ML stability prediction algorithms perform best should be answered decisively _before_ large DFT databases like MP or the OQMD commit significant resources to new efforts to expand their databases. In this work, we aim to answer which of these is the winning methodology in a future-proof benchmark that closely simulates using ML to guide a real-world discovery campaign.

We use the distance to the PBE DFT convex hull as the best proxy of true crystal stability available on large datasets, recognizing that our ground truth is ignorant of entropic stabilization or metastable states. Dealing with this significant layer of additional complexity is the premise of recent efforts to predict synthesis reaction pathways @mcdermott_graph-based_2021 @aykol_rational_2021 @wen_chemical_2023. These algorithms stand to benefit from more efficient estimates of the reaction energy barriers. In the past, such energy estimates often required expensive ab-initio simulation. With ML stability predictions gaining in fidelity even in OOD settings while remaining much cheaper than DFT as shown in this benchmark, we believe future efforts in reaction pathway finding can leverage the same ML models for similar acceleration factors.

## Related Work

### Using Stability rather than Formation Energies

In 2020, Chris Bartel et al. @bartel_critical_2020 benchmarked 7 models, finding all of them able to predict DFT formation energies with useful accuracy.
However, when asked to predict stability (specifically decomposition enthalpy), the performance of all models deteriorated sharply.
This insight meant that ML models are much less useful than DFT for discovering new solids than prior studies had suggested.
The paper identified two main reasons for the sharp decline in predictive power:

1. Stability is a property not only of the material itself but also the chemical space of competing phases it occupies. Current ML algorithms are given an input that only describes the single material they are asked to predict, leaving them clueless of competing phases.
1. Unlike DFT, ML models appear to benefit less from systematic error cancellation across similar chemistries. A first-principles theory of physics incurs similar errors for similar systems. In particular, the error direction (over- or underestimating the energy) should tend to agree across members of a chemical space. When looking at the relative energy differences that determine stability, systematic errors more strongly cancel. ML errors appear to follow a more random distribution, making cancellations less likely.

Bartel et al. stressed that to demonstrate the utility of ML for materials discovery, the vanity metric of formation energy accuracy must be replaced with stability predictions.
Moreover, the qualitative leap in performance from Roost @goodall_predicting_2020, the best compositional model benchmarked, to CGCNN @xie_crystal_2018, the single structural model they tested, shows structure plays a crucial role in determining the stability of materials.
However, using the DFT-relaxed structure as input to CGCNN renders the discovery pipeline circular as the input becomes the very thing we aim to find. A true test of prospective utility requires using unrelaxed structures as the next most information-rich input.

### Matbench

As the name suggests, this work seeks to expand upon the original Matbench suite of property prediction tasks @dunn_benchmarking_2020. By providing a standardized collection of datasets along with canonical cross-validation splits for model evaluation, Matbench helped focus the field of ML for materials, increase comparability across papers and provide a quantitative measure of progress in the field. It was a similar attempt to catalyze the field of ML for materials through competition and setting goal posts as ImageNet was for computer vision.

Matbench released a test suite of 13 supervised tasks for different material properties ranging from thermal (formation energy, phonon frequency peak), electronic (band gap), optical (refractive index) to tensile and elastic (bulk and shear moduli).
They range in size from ~300 to ~132,000 samples and include both DFT and experimental data sources. 4 tasks are composition-only while 9 provide the relaxed crystal structure as input.
Importantly, all tasks were exclusively concerned with the properties of known materials.
We believe a task that simulates a materials discovery campaign by requiring materials stability prediction from unrelaxed structures to be a missing piece here.

### The Open Catalyst Project

The Open Catalyst Project (OCP) is a large-scale initiative to discover substrate-adsorbate combinations that catalyze key industrial reactions which process said adsorbates into more useful products.
The OCP has released two data sets thus far, OCP20 @chanussot_open_2021 and OCP22 @tran_open_2022, for training and benchmarking ML models.

However, the ambition and scale of OCP comes with limitations for its use as a benchmark.
OCP20 is already 10x larger than the largest crystal structure data sets available, imposing a high barrier to entry for researchers without access to cloud-scale computing. Moreover, analyzing the learning behavior of the baseline models, the OCP authors estimate that up to 10 orders of magnitude more data (not 10x mind you) will be required before existing ML models reach the accuracy of DFT for adsorbate energies and become useful (see fig. 7 right @chanussot_open_2021).

In contrast, we believe the discovery of stable materials is a problem where ML methods have matured enough to be usefully deployed at scale after training existing models for only $\mathcal{O}(10^2)$ GPU hours.

## Data Sets

The choice of data for the train and test sets of this benchmark fell on the latest Materials Project (MP) @jain_commentary_2013 database release ([v2022.10.28](https://docs.materialsproject.org/changes/database-versions#v2022.10.28) at time of writing) and the WBM dataset @wang_predicting_2021.

### Materials Project Training Set

The Materials Project is a well-known effort to calculate the properties of all inorganic materials using high-throughput ab-initio methods. Seeded from a subset of the Inorganic Crystal Structure Database (ICSD) @allen_crystallographic_1999, the initial release of the database consisted of ~9 k crystals.
At time of writing, the Materials Project database has grown to [~154 k crystals](https://materialsproject.org/materials), covering diverse chemistries (see [periodic table heatmap](/about-the-data#--chemical-diversity)) and providing relaxed and initial structure as well as the relaxation trajectory for every entry.
For our benchmark, the training set is all data available from the 2022.10.28 MP release. Models are free to train on relaxed and/or unrelaxed structures or the full DFT relaxation trajectory. This flexibility is intended to allow authors to experiment and exploit the large variety of available data. It also aligns with our expectation that future progress in ML for crystal stability is more likely to result from better-leveraging training data than innovations in model architecture.

### WBM Test Set

The WBM data set @wang_predicting_2021 consists of ~257 k structures generated via chemical similarity-based elemental substitution of MP source structures followed by DFT relaxation and calculating each crystal's convex hull distance. Which element replaces an existing one in a given source structure was determined by randomly sampling according to the weights in a chemical similarity matrix mined from the ICSD. That is, elements are more likely to be replaced by elements that frequently co-occur in the ICSD for the given structure's prototype. But in principle an element could be replaced by any of the 89 elements present in MP (atomic numbers 1 - 84 and 89 - 94).

The WBM authors performed 5 iterations of this substitution process (we also refer to these as batches). After each step, the newly generated structures found to be stable after DFT relaxation flow back into the source pool to partake in the next round of substitution. This split of the data into batches of increasing substitution count is a unique and compelling feature of the test set as it allows out-of-distribution (OOD) testing by seeing if model performance degrades with substitution count. A higher number of elemental substitutions on average carries the structure further away from the region of material space covered by the MP training set. See [per-batch figures in the SI](/si#wbm-batch-robustness-as-a-measure-of-extrapolation-prowess) for details.

Throughout this work, we define stability as being on or below the convex hull of the MP training set. ~42 k out of ~257 k materials in WBM satisfy this criterion. Our code treats the stability threshold as a dynamic parameter for future more detailed performance analysis at different thresholds. For initial analysis in this direction, see ROC curves in [the SI](/si#roc-curves).

The WBM test set has an energy above the MP convex hull distribution with mean ± std = 0.02 ± 0.25 eV/atom (see [target distribution](https://matbench-discovery.janosh.dev/about-the-data#--target-distribution)). As WBM explores regions of materials space not well sampled by MP, many of the discovered materials that are stable w.r.t. MP's convex hull are not stable w.r.t. each other. Less than half or around ~20 k remain on the convex hull when merging the MP and WBM hulls, suggesting many WBM structures are repeated samples into the same new chemical spaces.
This observation highlights a critical aspect of this benchmark in that we purposely operate with an incomplete convex hull. Only current knowledge is accessible to a real discovery campaign and our metrics are designed to reflect this.

To simulate a real discovery campaign, our test set inputs are unrelaxed structures (obtained from element substitution on MP source structures, as mentioned above) while our target labels are the relaxed PBE formation energies. This task type was coined IS2RE (initial structure to relaxed energy) by OCP @chanussot_open_2021.

## Models

Our initial benchmark release includes 8 models. @Fig:model-metrics includes all models but we focus on the 6 best performers in subsequent figures for visual clarity.

1. **Voronoi+RF** @ward_including_2017 - A random forest trained to map a combination of composition-based Magpie features and structure-based relaxation-invariant Voronoi tessellation features (effective coordination numbers, structural heterogeneity, local environment properties, ...) to DFT formation energies.

   This old model predates most deep learning for materials but significantly improved over Coulomb matrix and partial radial distribution function methods. It therefore serves as a good baseline model to see how much value deep learning models are able to extract from the increasingly large training data on offer in this field.

1. **Wrenformer** - For this benchmark, we introduce Wrenformer which is a variation on Wren @goodall_rapid_2022 constructed using standard QKV-Transformer blocks @vaswani_attention_2017 in place of message-passing layers to reduce memory usage, allowing it to scale to structures with >16 Wyckoff positions.

   Like its predecessor, Wrenformer is a fast coordinate-free model aimed at accelerating screening campaigns where even the unrelaxed structure is a priori unknown.

   The key idea is that by training on the Wyckoff positions (symmetry-related positions in the crystal structure), the model learns to distinguish polymorphs while maintaining discrete and computationally enumerable inputs. During screening, we can take advantage of the fact that across the 230 crystallographic space groups in 3D, there are only 1731 different Wyckoff positions. This allows a fast model like Wrenformer to predict, for a given composition, the energy of all possible combinations of spacegroup and Wyckoff positions. A useful model should rank the symmetry exhibited by the actual crystal as one of its lowest predicted energies. Knowing what symmetry are present in a structure and how they confine atoms to planes, lines or even points then allows for a cheaper DFT relaxation with fewer degrees of freedom to obtain the relaxed structure.

1. **CGCNN** @xie_crystal_2018 - The Crystal Graph Convolutional Neural Network (CGCNN) was the first neural network model to directly learn 8 different DFT-computed material properties from a graph representing the atoms and bonds in a periodic crystal.

   CGCNN was among the first to show that just like in other areas of ML, given large enough training sets, neural networks can learn embeddings that reliably outperform all human-engineered structure features directly from the data.

1. **CGCNN+P** @gibson_data-augmentation_2022 - This work proposes a simple, physically motivated structure perturbations to augment stock CGCNN's training data of relaxed structures with structures resembling unrelaxed ones but mapped to the same DFT final energy. Here we chose $P=5$, meaning the training set was augmented with 5 random perturbations of each relaxed MP structure mapped to the same target energy.

   In contrast to all other structure-based GNNs considered in this benchmark, CGCNN+P is not attempting to learn the Born-Oppenheimer potential energy surface. The model is instead taught the PES as a step-function that maps each valley to its local minimum. The idea is that during testing on unrelaxed structures, the model will predict the energy of the nearest basin in the PES. The authors confirm this by demonstrating a lowering of the energy error on unrelaxed structures.

1. **MEGNet** @chen_graph_2019 - MatErials Graph Network is another GNN similar to CGCNN for material properties of relaxed structures that also updates the edge and global features (like pressure, temperature, entropy) in its message passing operation.

   This work showed that learned element embeddings encode periodic chemical trends and can be transfer-learned from large data sets (formation energies) to predictions on small data properties (band gaps, elastic moduli).

1. **M3GNet** @chen_universal_2022 - M3GNet is a GNN-based universal (as in full periodic table) interatomic potential (IAP) for materials trained on up to 3-body interactions in the initial, middle and final frame of MP DFT relaxations. The model takes the unrelaxed input and emulates structure relaxation before predicting energy for the pseudo-relaxed structure.

1. **M3GNet + MEGNet** @chen_universal_2022 @chen_graph_2019 - This combination of models uses M3GNet to relax initial structures and then passes it to MEGNet to predict the formation energy.

1. **BOSWR + MEGNet** @zuo_accelerating_2021 - BOWSR uses a symmetry-constrained Bayesian optimizer (BO) with a surrogate energy model (here MEGNet) to perform an iterative exploration-exploitation-based search of the potential energy landscape. The high sample count needed to explore the PES with BO makes this by far the most expensive model.

   The high sample count needed to explore the high-dimensional PES makes this type of "DFT-free" structure relaxation by far the most expensive model.

## Results

### Summary Table

<MetricsTable />

> @label:fig:model-metrics Regression and classification metrics for all models tested on our benchmark. The heat map ranges from yellow (best) to blue (worst) performance. DAF = discovery acceleration factor (see text), TPR = true positive rate, FNR = false negative rate, MAE = mean absolute error, RMSE = root mean squared error

@Fig:model-metrics shows performance metrics for all models considered in v1 of our benchmark.
M3GNet takes the top spot on most metrics and emerges as current SOTA for ML-guided materials discovery. The discovery acceleration factor (DAF) measures how many more stable structures a model found among the ones it predicted stable compared to the dummy discovery rate of 43 k / 257 k $\approx$ 16.7% achieved by randomly selecting test set crystals. The dummy MAE of always predicting the mean distance to the convex hull is 0.17 eV/atom. Our baseline model Voronoi RF beats random performance in both DAF and MAE.

The maximum possible DAF on our current test set is $\frac{1}{0.167} \approx 6$. This highlights the fact that our benchmark is made more challenging by deploying models on an already enriched space with a much higher fraction of stable structures compared to materials space at large. As the convex hull becomes more thoroughly sampled by future discovery, the fraction of unknown stable structures decreases, naturally leading to less enriched future test sets which will allow for higher maximum DAFs.

We note that while F1 score and DAF of models that make one-shot predictions directly from unrelaxed inputs (CGCNN, MEGNet, Wrenformer) are seemingly unaffected, the $R^2$ of these models is significantly worse. The reason MEGNet outperforms M3GNet on DAF becomes clear from @fig:cumulative-clf-metrics. MEGNet's line ends at 41.6 k materials which is closest to the true number of 43 k stable materials. All other models overpredict this number by anywhere from 40% (~59 k for CGCNN) to 104% (87 k for Wrenformer), resulting in large numbers of false positive predictions that drag down their DAFs.

### Cumulative Classification Metrics

{#if browser}
<CumulativeClfMetrics style="margin: 0 -2em 0 -4em;" />
{/if}

> @label:fig:cumulative-clf-metrics Running precision and recall over the course of a simulated discovery campaign. This figure highlights how different models perform better or worse depending on the length of the discovery campaign. Length here is an integer measuring how many DFT relaxations you have compute budget for.

@Fig:cumulative-clf-metrics simulates ranking materials from most to least stable according to model predictions and going down the list calculating the precision and recall of correctly identified stable materials at each step, i.e. exactly how these models could be used in a prospective materials discovery campaign.

A line terminates when a model believes there are no more materials in the WBM test set below the MP convex hull. The dashed vertical line shows the actual number of materials on or below the MP hull. Most models overestimate the number of stable materials. The dashed diagonal Optimal Recall line would be achieved if a model never made a false negative prediction and predicts everything as unstable exactly when the true number of stable materials is reached. Zooming in on the top-left corner of the precision plot, we observe that MEGNet and the combo M3GNet + MEGNet are particularly suitable for very short discovery campaigns of less than 2000 and 3000 structure relaxations, respectively, after which M3GNet takes the lead on both precision and recall and does not relinquish for the remainder.

### Rolling MAE vs. Hull Distance

{#if browser}
<RollingMaeModels />
{/if}

> @label:fig:rolling-mae-vs-hull-dist-models Rolling MAE on the WBM test set as the energy to the convex hull of the MP training set is varied. The width of the box in the bottom corner indicates the rolling window within which errors were averaged (think smoothing strength). The red-highlighted 'triangle of peril' shows where the models are most likely to misclassify structures. As long as a model's rolling MAE remains inside the triangle, its mean error is larger than the distance to the convex hull. If the model's error happens to point towards the stability threshold at 0 eV from the hull (the plot's center), it's average error will change the stability classification of a material from true positive/negative to false negative/positive.

@Fig:rolling-mae-vs-hull-dist-models visualizes a model's reliability as a function of a material's hull distance. The lower its rolling MAE exits the shaded triangle, the better. Inside this area, the model's mean error is larger than the distance to the convex hull, making misclassifications likely. Outside the triangle even if the model's error points toward the stability threshold at 0 eV from the hull (the plot's center), the mean error is too small to move a material over the stability threshold which would cause a false stability classification. M3GNet achieves the lowest overall MAE and exits the peril zone much sooner than other models on the right half of the plot. This means it rarely misclassifies unstable materials that lie more than 40 meV above the hull. On the plot's left half, CGCNN+P exits the peril zone first, albeit much further from the hull at more than 100 meV below. Essentially, all models are prone to false negative predictions even for materials far below the known hull which aligns with the smaller amount of training data for materials on or below the known convex hull (see the [test set's target distribution](/about-the-data#--target-distribution)).

### Predicted Hull Distance Parity Plots

{#if browser}
<EachScatter />
{/if}

> @label:fig:each-scatter-models Parity plot for each model's energy above hull predictions (based on their formation energy preds) vs DFT ground truth

TODO: mention we consistently see deducting old MP corrections and applying 2020 scheme from MEGNEt e_form predictions increases MAE, no matter if paired with BOWSR, M3GNet, CHGNet or standalone

### Per-Element Model Error Heatmap

<ElementErrorsPtableHeatmap current_model={["Mean over models"]} />

We observe that many models - notably excluding the ML-IAPs CHGNet and M3GNet - struggle with the halogens. The fact that this increase in error is largest for fluorine, the halogen with the largest electronegativity, suggests a relation to the high degree of ionic bonds in materials containing halogens. GNNs may face difficulties in accurately capturing the long-range interactions required to describe such ionic bonds. However, Wrenformer and Voronoi RF do not include bond lengths in their input features yet also exhibit increased errors on halogens.

We don't believe this is a simple case of wider energy ranges for halogens in the test set making the prediction task inherently more difficult. Halogen errors remain elevated even when normalizing for data range by dividing each element's error by the standard deviation over target energies of all structures containing said element.

Many GNN message-passing functions incorporate a soft attention coefficient designed to evaluate the significance of one atom's contribution to another. Insufficient training support for materials with ionic bonds could result in this coefficient primarily learning covalent interactions. Considering fluorine is the 10th most abundant element in MP (~12k structures) with other halogens not far behind, this also seems unlikely. As such, we believe this phenomenon invites further investigation into possible pitfalls in ionic material stability prediction.

## Discussion

From @fig:model-metrics we see several models achieve a DAF > 2 in this realistic benchmark scenario.
Consequently, the benefits of deploying ML-based triage in high-throughput computational materials discovery applications likely warrant the time and setup required.
However, there are many aspects on which further progress is necessary, for example, models still make large numbers of false positive predictions for materials over 50 meV above the convex hull and much less likely to be synthesizable, greatly reducing the DAF.
The results obtained from version 1 of our benchmark show that ML universal interatomic potentials like M3GNet are the most promising methodology to pursue going forward, being both ~20x cheaper to run than black box optimizers like BOWSR and having access to more training structures than coordinate-free approaches like Wrenformer.

<!-- Before obtaining the results for this benchmark, we saw grounds for debate on which is the most promising ML discovery methodology to develop further. After all, some authors of this work were involved with developing Wren/Wrenformer, a one-shot energy predictor using only coarse-grained symmetry features to predict final energy from initial structure. -->

Just like minerals, oil or any other finite resource, the task of discovering new materials will necessarily become more challenging over time as currently undersampled regions of materials space are explored. Before the results of this benchmark, we saw grounds for debate on which is the most promising ML discovery methodology to develop further. We now believe the path to making ML a ubiquitous discovery tool to be straightforward and one the field is already pursuing: training foundational IAPs on significantly more data may get us there even without further algorithmic or model improvements. If we manage to double or triple the discovery acceleration factor, the benefits of including ML in materials discovery workflows will significantly outweigh the setup costs.

We welcome further model submissions as well as data contributions for version 2 of this benchmark to the GitHub repo at
[{repo}]({repo}).

## Acknowledgments

Janosh Riebesell acknowledges support from the German Academic Scholarship Foundation ([Studienstiftung](https://wikipedia.org/wiki/Studienstiftung)) and gracious hosting as a visiting affiliate in the groups of Kristin Persson and Anubhav Jain.

We would like to thank Jason Gibson, Ekin Dogus Cubuk, Tian Xie, Chi Chen and Ryota Tomioka for helpful discussions.

## Author Contributions

Janosh Riebesell: Methodology, Software, Data Curation, Formal analysis. Rhys Goodall: Conceptualization, Software, Formal analysis. Anubhav Jain: Supervision. Kristin Persson: Supervision. Alpha Lee: Supervision.
