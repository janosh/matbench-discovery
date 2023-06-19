<script>
  import { repository as repo } from '$site/package.json'
  import EachScatterModels from '$figs/each-scatter-models-3x3.svelte'
  import MetricsTable from '$figs/metrics-table.svelte'
  import CumulativePrecisionRecall from '$figs/cumulative-precision-recall.svelte'
  import RollingMaeVsHullDistModels from '$figs/rolling-mae-vs-hull-dist-models.svelte'
  import ElementErrorsPtableHeatmap from '$models/element-errors-ptable-heatmap.svelte'
  import HistClfTrueHullDistModels from '$figs/hist-clf-pred-hull-dist-models-3x3.svelte'
  import { onMount } from 'svelte'

  let mounted = false
  onMount(() => (mounted = true))
</script>

<summary class="abstract">

We present a new machine learning (ML) evaluation framework for materials stability predictions named `Matbench Discovery`. Our task closely simulates the deployment of ML energy models in a high-throughput search for stable inorganic crystals. It is accompanied by an interactive leaderboard and a Python package for easy ingestion of our training/test sets into future model submissions. To answer the question which ML methodology performs best at materials discovery, we explore a wide variety of models. Our initial selection ranges from random forests to GNNs, from one-shot predictors to iterative Bayesian optimizers and universal interatomic potentials (UIP) that closely emulate DFT. We find UIPs to be in a class of their own, achieving the highest F1 scores and discovery acceleration factors (DAF) of more than 3, i.e. 3x more stable structures found compared to dummy selection in our already enriched search space. We also identify a sharp disconnect between commonly used regression metrics and more task-relevant classification metrics. CGCNN and MEGNet are worse than dummy regressors, but substantially better than dummy classifiers, suggesting that the field overemphasizes the wrong performance indicators. Our results highlight the need to optimize metrics that measure true stability hit rate improvements and provide valuable insights for maintainers of high throughput materials databases by demonstrating that these models have matured enough to play a vital role as pre-filtering steps to effectively allocate compute budget for DFT relaxations.

</summary>

## Introduction

### Why Materials Discovery?

@davies_computational_2016 identified $~10^{10}$ possible quaternary materials allowed by electronegativity and charge-balancing rules. Of these, only $~10^5$ are known experimentally and $~10^6$ have been simulated. The space of quinternaries and higher is even less explored, leaving vast numbers of potentially useful materials to be discovered. The discovery of new materials is a key driver of technological progress and lies on the path to more efficient solar cells, lighter and longer-lived batteries, smaller and more efficient transistor gates just to name a few. The materials we have mastered define our technological abilities so much so that we name eras of human development after them, going from stone, bronze and iron age to our current silicon age. In light of global warming, these advances cannot come fast enough. Any speed-up new methods might yield should be leveraged to their fullest extent.

### Why use machine learning for discovery?

Despite significant advances in empirical, theoretical and computational materials science, discovering new materials still requires labor-intensive experimentation, complex calculations, trial-and-error and often happens fortuitously rather than through rational design. ML methods can be useful for materials discovery because they efficiently extract and distill trends from huge datasets (think single GPTs able to encode all human languages), can handle high dimensionality, multiple objectives, uncertainty, and noisy or sparse data.

Statistical models are also orders of magnitude faster than ab-initio simulation. As such, they are most suitable for high-throughput searches as a pre-filtering step to more expensive, higher-fidelity simulation methods. The use of neural networks for learning the density-functional theory (DFT) potential energy surface (PES) can be traced as far back as @behler_generalized_2007. This work opened the floodgates for material scientists to devote significant effort to fitting ever more sophisticated models to known samples of the PES.
Initially, most of these models were trained and deployed as interatomic potentials to study known materials of interest which required curating custom training data time and time again for each application @bartok_machine_2018 @deringer_general-purpose_2020.
As larger and more diverse datasets emerged from initiatives like the Materials Project (MP) @jain_commentary_2013 or the Open Quantum Materials Database (OQMD) @saal_materials_2013, researchers have begun to train models that cover the full periodic table, opening up the prospect of ML-guided materials discovery to increase hit rate and speed of DFT/expert-driven searches.

### Limitations of current benchmarking in ML

Yet despite many advances in ML for materials, it is unclear which methodology performs best at predicting material stability, let alone which model. The field lacks a standardized benchmark task that accurately simulates applying models in a prospective materials discovery campaign.

1. **Lack of realism**: Benchmark tasks can be idealized and simulate overly simplified conditions that do not reflect the real-world challenges a model faces when used in an actual discovery campaign. This can lead to pretty leaderboards listing SOTA models that underwhelm when used in production and has caused some disillusionment with ML methods in the past. Examples of how this comes about are choosing the wrong target or picking an unrepresentative train/test split.
1. **Limited diversity**: Benchmark datasets may be too small and contain only a limited number of materials, unrepresentative of the huge diversity of materials space. This can make models look good even if they fail to generalize.
1. **Opportunity cost**: Some benchmarks give insufficient thought to the high cost of failed experiments. Looking purely at global metrics like $\text{MAE}$, $\text{RMSE}$ and $R^2$ can give practitioners a false sense of security. Accurate regressors are susceptible to unexpectedly high false-positive rates if those accurate predictions lie close to the decision boundary. This can cause experimentalists to waste time and resources on doomed experiments.
<!-- TODO who to cite here -->
1. **Scalability**: Many benchmark tasks have too small data to adequately simulate the high-throughput and large-date regimes that future discovery efforts are likely to encounter. Confining model testing to the small data regime can obfuscate poor scaling relations like Gaussian Processes whose training costs grow cubically with training sample count or random forests that achieve outstanding performance on few data points but fail to extract the full information content out of larger datasets, leading to flatter learning curves compared to neural networks and worse performance when large amounts of training data are available.

Ideally, the question of which ML stability prediction algorithms perform best should be answered decisively _before_ large DFT databases like MP or the OQMD commit resources to ML efforts to expand their databases. In this work, we aim to answer which of these is the winning methodology in a future-proof benchmark that closely simulates using ML to guide a real-world discovery campaign.

We use the distance to the PBE DFT convex hull as the best proxy of true crystal stability available on large datasets, recognizing that our ground truth is ignorant of entropic stabilization or metastable states. Dealing with this significant layer of additional complexity is the premise of recent efforts to predict synthesis reaction pathways @mcdermott_graph-based_2021 @aykol_rational_2021 @wen_chemical_2023. These algorithms stand to benefit from more efficient estimates of the reaction energy barriers. In the past, such energy estimates often required expensive ab-initio simulation. With ML stability predictions gaining in fidelity even on OOD samples while remaining much cheaper than DFT as shown in this benchmark, we believe future efforts in reaction pathway finding can leverage the same ML models for similar acceleration factors.

## Related Work

### Using Stability rather than Formation Energies

In 2020, Chris Bartel et al. @bartel_critical_2020 benchmarked 7 models, finding all of them able to predict DFT formation energies with useful accuracy.
However, when asked to predict stability (specifically decomposition enthalpy), the performance of all models deteriorated sharply.
This insight meant that ML models are much less useful than DFT for discovering new solids than prior studies had suggested.
The paper identified two main reasons for the sharp decline in predictive power:

1. Stability is a property not only of the material itself but also the chemical space of competing phases it occupies. Current ML algorithms are given an input that only describes the single material they are asked to predict, leaving models clueless of competing phases.
1. Unlike DFT, ML models appear to benefit less from systematic error cancellation across similar chemistries. A first-principles theory of physics incurs similar errors for similar systems. In particular, the error direction (over- or underestimating the energy) tends to agree across members of the same chemical space. When looking at the relative energy differences that determine stability, systematic errors more strongly cancel. ML errors appear to follow a more random distribution, making cancellations less likely.

Bartel et al. stressed that to demonstrate the utility of ML for materials discovery, the vanity metric of formation energy accuracy must be replaced with stability predictions.
Moreover, the qualitative leap in performance from Roost @goodall_predicting_2020, the best compositional model benchmarked, to CGCNN @xie_crystal_2018, the single structural model they tested, highlighted that structure plays a crucial role in determining the stability of materials.
However, using the DFT-relaxed structure as input to CGCNN renders the discovery pipeline circular as the input becomes the very thing we aim to find. A true test of prospective utility requires using unrelaxed structures as the next most information-rich input.

### Matbench

As the name suggests, this work seeks to expand upon the original Matbench suite of property prediction tasks @dunn_benchmarking_2020. By providing a standardized collection of datasets along with canonical cross-validation splits for model evaluation, Matbench helped focus the field of ML for materials, increase comparability across papers and provide a quantitative measure of progress in the field. It aimed to catalyze the field of ML for materials through competition and establishing common goal posts in a similar fashion as ImageNet did for computer vision.

Matbench released a test suite of 13 supervised tasks for different material properties ranging from thermal (formation energy, phonon frequency peak), electronic (band gap), optical (refractive index) to tensile and elastic (bulk and shear moduli).
They range in size from ~300 to ~132,000 samples and include both DFT and experimental data sources.
Importantly, all tasks were exclusively concerned with the properties of known materials.
We believe a task that simulates a materials discovery campaign by requiring materials stability prediction from unrelaxed structures to be a missing piece here.

### The Open Catalyst Project

The Open Catalyst Project (OCP) is a large-scale initiative to discover substrate-adsorbate combinations that catalyze key industrial reactions which process said adsorbates into more useful products.
The OCP has released two data sets thus far, OCP20 @chanussot_open_2021 and OCP22 @tran_open_2022, for training and benchmarking ML models.

However, the ambition and scale of OCP comes with limitations for its use as a benchmark.
OCP20 is already 10x larger than the largest crystal structure data sets available, imposing a high barrier to entry for researchers without access to cloud-scale computing. Moreover, analyzing the learning behavior of the baseline models, the OCP authors estimate that up to 10 orders of magnitude more data (not 10x, mind you) will be required before existing ML models reach the accuracy of DFT for adsorbate energies and become useful (see fig. 7 right @chanussot_open_2021).

In contrast, we believe the discovery of stable materials is a problem where ML methods have matured enough to be usefully deployed at scale after training existing models for only $\mathcal{O}(10^2)$ GPU hours.

## Data Sets

The choice of data for the train and test sets of this benchmark fell on the latest Materials Project (MP) @jain_commentary_2013 database release ([v2022.10.28](https://docs.materialsproject.org/changes/database-versions#v2022.10.28) at time of writing) and the WBM dataset @wang_predicting_2021.

### Materials Project Training Set

The Materials Project is a well-known effort to calculate the properties of all inorganic materials using high-throughput ab-initio methods. Seeded from a subset of the Inorganic Crystal Structure Database (ICSD) @allen_crystallographic_1999, the initial release of the database consisted of ~9 k crystals.
At time of writing, the Materials Project database has grown to [~154 k crystals](https://materialsproject.org/materials), covering diverse chemistries (see [periodic table heatmap](/about-the-data#--chemical-diversity)) and providing relaxed and initial structure as well as the relaxation trajectory for every entry.

Our benchmark defines the training set as all data available from the [2022.10.28 MP release](https://docs.materialsproject.org/changes/database-versions#v2022.10.28). Models are free to train on relaxed and/or unrelaxed structures in MP or the full DFT relaxation trajectory. We recorded a snapshot of energies, forces, stresses and magnetic moments for all MP ionic steps on 2023-03-15 as the canonical training set for v1 of Matbench Discovery, and provide convenience functions through our [Python package](https://pypi.org/project/matbench-discovery) for easily feeding that data into future model submissions to our benchmark. This flexibility is intended to allow authors to experiment with and exploit the large variety of available data. It also aligns with our expectation that future progress in ML for crystal stability is more likely to result from better leveraging training data than innovations in model architecture.

### WBM Test Set

The WBM data set @wang_predicting_2021 consists of ~257 k structures generated via chemical similarity-based elemental substitution of MP source structures followed by DFT relaxation and calculating each crystal's convex hull distance. Which element replaces an existing one in a given source structure was determined by randomly sampling according to the weights in a chemical similarity matrix mined from the ICSD. That is, elements are more likely to be replaced by elements that frequently co-occur in the ICSD for the given structure prototype. But in principle an element could be replaced by any of the 89 elements present in MP (atomic numbers 1 - 84 and 89 - 94).

The WBM authors performed 5 iterations of this substitution process (we refer to these as batches). After each step, the newly generated structures found to be thermodynamically stable after DFT relaxation flow back into the source pool to partake in the next round of substitution. This split of the data into batches of increasing substitution count is a unique and compelling feature of the test set as it allows out-of-distribution (OOD) testing by seeing if model performance degrades with substitution count. A higher number of elemental substitutions on average carries the structure further away from the region of material space covered by the MP training set. See [per-batch figures in the SI](/si#wbm-batch-robustness-as-a-measure-of-extrapolation-prowess) for details.

Throughout this work, we define stability as being on or below the convex hull of the MP training set. ~42 k out of ~257 k materials in WBM satisfy this criterion. Our code treats the stability threshold as a dynamic parameter for future more detailed performance analysis at different thresholds. For initial analysis in this direction, see ROC curves in [the SI](/si#roc-curves).

The WBM test set has an energy above the MP convex hull distribution with mean ± std = 0.02 ± 0.25 eV/atom (see [target distribution](https://janosh.github.io/matbench-discovery/about-the-data#--target-distribution)). As WBM explores regions of materials space not well sampled by MP, many of the discovered materials that are stable w.r.t. MP's convex hull are not stable w.r.t. each other. Less than half or around ~20 k remain on the convex hull when merging the MP and WBM hulls, suggesting many WBM structures are repeated samples into the same new chemical spaces.
This observation highlights a critical aspect of this benchmark in that we purposely operate with an incomplete convex hull. Only current knowledge is accessible to a real discovery campaign and our metrics are designed to reflect this.

To simulate a real discovery campaign, our test set inputs are unrelaxed structures (obtained from element substitution on MP source structures, as mentioned above) while our target labels are the relaxed PBE formation energies. This task type was coined IS2RE (initial structure to relaxed energy) by OCP @chanussot_open_2021.

## Models

Our initial benchmark release includes 8 models.

1. **Voronoi+RF** @ward_including_2017 - A random forest trained to map a combination of composition-based Magpie features and structure-based relaxation-invariant Voronoi tessellation features (effective coordination numbers, structural heterogeneity, local environment properties, ...) to DFT formation energies.

   This older model predates most deep learning for materials but significantly improved over Coulomb matrix and partial radial distribution function methods. It therefore serves as a good baseline model to see how much value deep learning models are able to extract from the increasingly large training data on offer in this field.

1. **Wrenformer** - For this benchmark, we introduce Wrenformer which is a variation on Wren @goodall_rapid_2022 constructed using standard QKV-Transformer blocks @vaswani_attention_2017 in place of message-passing layers to reduce memory usage, allowing it to scale to structures with >16 Wyckoff positions.

   Like its predecessor, Wrenformer is a fast coordinate-free model aimed at accelerating screening campaigns where even the unrelaxed structure is a priori unknown.

   The key idea is that by training on the Wyckoff positions (symmetry-related positions in the crystal structure), the model learns to distinguish polymorphs while maintaining discrete and computationally enumerable inputs. During screening, we can take advantage of the fact that across the 230 crystallographic space groups in 3D, there are only 1731 different Wyckoff positions. This allows a fast model like Wrenformer to predict, for a given composition, the energy of all possible combinations of spacegroup and Wyckoff positions for a maximum unit cell size allowing for a single symmetry-constrained relaxation of the lowest ranked prototype to be performed or selection of an enriched subset of prototypes as the initial population for a crystal structure searching algorithm.

1. **CGCNN** @xie_crystal_2018 - The Crystal Graph Convolutional Neural Network (CGCNN) was the first neural network model to directly learn 8 different DFT-computed material properties from a graph representing the atoms and bonds in a periodic crystal.

   CGCNN was among the first to show that just like in other areas of ML, given large enough training sets, neural networks can learn embeddings that reliably outperform all human-engineered structure features directly from the data.

1. **CGCNN+P** @gibson_data-augmentation_2022 - Identical to CGCNN except for a difference in training procedure. The +P stands for training set augmentation using random structure perturbations. We apply a slight lattice strain and nudge the atoms but keep the same energy target to create additional training samples. Here we chose $P=5$ meaning random perturbations are repeated 5 times for each relaxed MP structure.

   In contrast to all other structure-based GNNs considered in this benchmark, CGCNN+P is not attempting to learn the Born-Oppenheimer potential energy surface. The model is instead taught the PES as a step-function that maps each valley to its local minimum. The idea is that during testing on unrelaxed structures, the model will predict the energy of the nearest basin in the PES. The authors confirm this by demonstrating a lowering of the energy error on unrelaxed structures.

1. **MEGNet** @chen_graph_2019 - MatErials Graph Network is another GNN similar to CGCNN for material properties of relaxed structures that also updates the edge and global features (like pressure, temperature, entropy) in its message passing operation.

   This work showed that learned element embeddings encode periodic chemical trends and can be transfer-learned from large data sets (formation energies) to predictions on small data properties (band gaps, elastic moduli).

1. **M3GNet** @chen_universal_2022 - M3GNet is a GNN-based universal (as in full periodic table) interatomic potential (UIP) for materials trained on up to 3-body interactions in the initial, middle and final frame of MP DFT relaxations. The model takes the unrelaxed input and emulates structure relaxation before predicting energy for the pseudo-relaxed structure.

1. **CHGNet** @deng_chgnet_2023 - CHGNet is another UIP for charge-informed atomistic modeling. It's distinguishing feature is that it was trained to predict magnetic moments on top of energies, forces and stresses in the MPTraj dataset consisting of relaxation trajectories for ~1.5 million MP structures. By modeling magnetic moments, CHGNet learns to accurately represent the orbital occupancy of electrons which allows it to predict both atomic and electronic degrees of freedom.

1. **BOSWR + MEGNet** @zuo_accelerating_2021 - BOWSR uses a symmetry-constrained Bayesian optimizer (BO) with a surrogate energy model (here MEGNet) to perform an iterative exploration-exploitation-based search of the potential energy landscape. The high sample count needed to explore the PES with BO makes this by far the most expensive model.

   The high sample count needed to explore the high-dimensional PES makes this type of "DFT-free" structure relaxation by far the most expensive model.

## Results

### Metrics Table

<MetricsTable />

> @label:fig:metrics-table Classification and regression metrics for all models tested on our benchmark. The heat map ranges from yellow (best) to blue (worst) performance. The dummy classifier uses the `scikit-learn` `stratified` strategy of randomly assigning stable/unstable labels according to the training set prevalence. The dummy regression metrics MAE, RMSE and $R^2$ are attained by always predicting the test set mean. Note that Voronoi RF, CGCNN and MEGNet are worse than dummy on regression metrics but better on some of the classification metrics, highlighting the importance of looking at the right metrics for the task at hand to gauge model performance.
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
CHGNet takes the top spot on all metrics except true positive rate (TPR) and emerges as current SOTA for ML-guided materials discovery. The discovery acceleration factor (DAF) measures how many more stable structures a model found compared to the dummy discovery rate of 43k / 257k $\approx$ 16.7\% achieved by randomly selecting test set crystals. Consequently, the maximum possible DAF is ~6. This highlights the fact that our benchmark is made more challenging by deploying models on an already enriched space with a much higher fraction of stable structures than uncharted materials space at large. As the convex hull becomes more thoroughly sampled by future discovery, the fraction of unknown stable structures decreases, naturally leading to less enriched future test sets which will allow for higher maximum DAFs.

Note that MEGNet outperforms M3GNet on DAF (2.70 vs 2.66) even though M3GNet is superior to MEGNet in all other metrics. The reason is the one outlined in the previous paragraph as becomes clear from @fig:cumulative-precision-recall. MEGNet's line ends at 55.6 k materials which is closest to the true number of 43 k stable materials in our test set. All other models overpredict the sum total of stable materials by anywhere from 40% (~59 k for CGCNN) to 104% (85 k for Wrenformer), resulting in large numbers of false positive predictions which lower their DAFs.

As noted, this is only a problem in practice for exhaustive discovery campaigns that validate _all_ stable predictions from a model. More frequently, model predictions will be ranked most-to-least stable and validation stops after some pre-determined compute budget is spent, say, 10k DFT relaxations. In that case, most of the false positive predictions near the less stable end of the candidate list are ignored and don't harm the campaign's overall discovery count.

We find a large performance gap between models that make one-shot predictions directly from unrelaxed inputs such as MEGNet, Wrenformer, CGCNN, CGCNN+P, Voronoi RF versus UIPs that predict forces to emulate DFT relaxation. While the F1 scores and DAFs of non-UIPs are seemingly unaffected, their $R^2$ coefficients are significantly worse. Except for CGCNN+P, all fail to achieve positive $R^2$. This means their predictions explain the observed variation in the data less than a horizontal line through the test set mean. In other words, these models are not predictive in a global sense (across the full dataset range). However, even models with negative $R^2$ can be locally good in the positive and negative tails of the test set hull distance distribution. They suffer most in the mode near the stability threshold of 0 eV/atom above the hull. This reveals an important shortcoming of $R^2$ as a metric for classification tasks like ours.

The reason CGCNN+P achieves better regression metrics than CGCNN but is still worse as a classifier becomes apparent from [the SI histograms](/si#fig:hist-clf-pred-hull-dist-models) by noting that the CGCNN+P histogram is more sharply peaked at the 0 hull distance stability threshold. This causes even small errors in the predicted convex hull distance to be large enough to invert a classification. Again, this is evidence to choose carefully which metrics to optimize. Regression metrics are far more prevalent when evaluating energy predictions. In our benchmark, energies are just means to an end to classify compound stability. Regression accuracy is of little use on its own unless it helps classification. The field needs to be aware that this is not a given.

### Cumulative Precision + Recall

{#if mounted}
<CumulativePrecisionRecall style="margin: 0 -2em 0 -4em;" />
{/if}

> @label:fig:cumulative-precision-recall Cumulative precision and recall over the course of a simulated discovery campaign. This figure highlights how different models perform better or worse depending on the length of the discovery campaign. Length here is an integer measuring how many DFT relaxations you have compute budget for. We only show the 6 best performing models for visual clarity.

@Fig:cumulative-precision-recall simulates ranking materials from most to least stable according to model-predicted energies. For each model, we go down that list material by material, calculating at each step the precision and recall of correctly identified stable materials. This simulates exactly how these models might be used in a prospective materials discovery campaign and reveal how a model's performance changes as a function of the discovery campaign length, i.e. the amount of resources available to validate model predictions.

A line terminates when a model believes there are no more materials in the WBM test set below the MP convex hull. The dashed vertical line shows the actual number of stable structures in our test set. All models are biased towards stability to some degree as they all overestimate this number, most of all Wrenformer by over 112%, least of all MEGNet by 30%. The diagonal Optimal Recall line would be achieved if a model never made a false negative prediction and stops predicting stable crystals exactly when the true number of stable materials is reached. Zooming in on the top-left corner of the precision plot, we observe that CHGNet is the only model without a sudden drop in precision right at the start of the discovery campaign. It keeps a strong lead over the runner-up M3GNet until reaching ~3k screened materials. From there, the CHGNet and M3GNet lines slowly converge until they almost tie at a precision of 0.52 after ~56k screened materials. At that point, CHGNet's list of stable predictions is exhausted while M3GNet continues to drop to 0.45 at 76 k, attributable to many false positives towards the end of the stable list.

### Rolling MAE vs. Hull Distance

{#if mounted}
<RollingMaeVsHullDistModels />
{/if}

> @label:fig:rolling-mae-vs-hull-dist-models Rolling MAE on the WBM test set as the energy to the convex hull of the MP training set is varied. The width of the box in the bottom corner indicates the rolling window within which errors were averaged (think smoothing strength). The red-highlighted 'triangle of peril' shows where the models are most likely to misclassify structures. As long as a model's rolling MAE remains inside the triangle, its mean error is larger than the distance to the convex hull. If the model's error happens to point towards the stability threshold at 0 eV from the hull (the plot's center), it's average error will change the stability classification of a material from true positive/negative to false negative/positive.

@Fig:rolling-mae-vs-hull-dist-models visualizes a model's reliability as a function of a material's hull distance. The lower its rolling MAE exits the shaded triangle, the better. Inside this area, the model's mean error is larger than the distance to the convex hull, making misclassifications likely. Outside the triangle even if the model's error points toward the stability threshold at 0 eV from the hull (the plot's center), the mean error is too small to move a material over the stability threshold which would cause a false stability classification. M3GNet achieves the lowest overall MAE and exits the peril zone much sooner than other models on the right half of the plot. This means it rarely misclassifies unstable materials that lie more than 40 meV above the hull. On the plot's left half, CGCNN+P exits the peril zone first, albeit much further from the hull at more than 100 meV below. Essentially, all models are prone to false negative predictions even for materials far below the known hull which aligns with the smaller amount of training data for materials on or below the known convex hull (see the [test set's target distribution](/about-the-data#--target-distribution)).

### Classification Histograms

{#if mounted}
<HistClfTrueHullDistModels />
{/if}

> @label:fig:hist-clf-true-hull-dist-models These histograms show the classification performance of models as a function of model-predicted hull distance on the $x$ axis. Models are sorted top to bottom by F1 score. While CHGNet and M3GNet perform almost equally well overall, these plots reveal that they do so via different trade-offs. M3GNet commits fewer false negative but more false positives predictions compared to CHGNet. In a real discovery campaign, false positives have a higher opportunity cost than false negatives since they result in wasted DFT relaxations or even synthesis time in the lab. A false negative by contrast is just one missed opportunity out of many.
> This observation is also reflected in the higher TPR and lower TNR of M3GNet vs CHGNet in @fig:metrics-table, as well as the lower error for CHGNet vs M3GNet on the stable side (left half) of @fig:rolling-mae-vs-hull-dist-models and M3GNet over CHGNet on the unstable side (right half) of @fig:rolling-mae-vs-hull-dist-models.

### Predicted Hull Distance Parity Plots

{#if mounted}
<EachScatterModels />
{/if}

> @label:fig:each-scatter-models Parity plot of DFT hull distance vs model hull distance predictions (derived from predicted formation energies). All models do well on the outliers. They suffer most in the mode of the distribution around the convex hull.
> Interestingly, all models do very well on the outliers. Where they suffer is in the mode of the distribution near the convex hull. All models exhibit a horizontal spike at 0 predicted hull distance for crystals that are actually very unstable, resulting in false positive predictions. Some models, Wrenformer in particular, also have a spike pointing upwards, which are materials actually close to the hull but predicted to be highly unstable by the model. These are false negatives, or missed needles in the haystack we're searching that is materials space.

<!-- TODO maybe mention we consistently see deducting old MP corrections and applying 2020 scheme from MEGNEt e_form predictions increases MAE, no matter if paired with BOWSR, M3GNet, CHGNet or standalone -->

### Per-Element Model Error Heatmap

<ElementErrorsPtableHeatmap current_model={["Mean error all models"]} />

We observe that many models - notably excluding the UIPs CHGNet and M3GNet - struggle with the halogens. The fact that this increase in error is largest for fluorine, the halogen with the largest electronegativity, suggests a relation to the high degree of ionic bonds in materials containing halogens. GNNs may face difficulties in accurately capturing the long-range interactions required to describe such ionic bonds. However, Wrenformer and Voronoi RF do not include bond lengths in their input features yet also exhibit increased errors on halogens.

We don't believe this is a simple case of wider energy ranges for halogens in the test set making the prediction task inherently more difficult. Halogen errors remain elevated even when normalizing for data range by dividing each element's error by the standard deviation over target energies of all structures containing said element.

Many GNN message-passing functions incorporate a soft attention coefficient designed to evaluate the significance of one atom's contribution to another. Insufficient training support for materials with ionic bonds could result in this coefficient primarily learning covalent interactions. Considering fluorine is the 10th most abundant element in MP (~12k structures) with other halogens not far behind, this also seems unlikely. As such, we believe this phenomenon invites further investigation into possible pitfalls in ionic material stability prediction.

<!-- TODO consider reporting that Wrenformer has higher batch variability in its rolling MAE than MEGNet. maybe sth to say for relaxation enabling more robust extrapolation -->

## Discussion

@Fig:metrics-table shows that several models achieve a DAF (Discovery Acceptance Factor) greater than 2.5 in a realistic benchmark scenario. CHGNet, which is currently the best-performing model, achieves a DAF of 3.06. These results demonstrate the effectiveness of using machine learning-based triage in high-throughput computational materials discovery applications. We contend our results demonstrate the real-world utility UIPs like CHGNet and M3GNet have attained for materials discovery, perhaps without yet having entered into common knowledge. We therefore recommend all future discovery efforts invest the time and resources necessary to leverage these new tools to their full potential.
However, many aspects remain on which further progress is possible. For example, even the best models in our analysis still make large numbers of false positive predictions for materials even 50 meV or more above the convex hull. These materials are thermodynamically unstable enough to make successful synthesis unlikely. When taking a prediction into the lab, this would negatively affect the real-world discovery acceleration factor.
The results obtained from version 1 of our benchmark show that ML universal interatomic potentials like M3GNet are the most promising methodology to pursue going forward, being both $\sim$20x cheaper to run than black box optimizers like BOWSR while also able to extract a more accurate and generalizable map of the PES given the same number of training structures than coordinate-free approaches like Wrenformer.

Just like minerals, oil or any other finite resource, the task of discovering new materials will necessarily become more challenging over time as currently undersampled regions of materials space are explored. Before the results of this benchmark, we saw grounds for debate on which is the most promising ML discovery methodology to develop further. We now believe the path to making ML a ubiquitous discovery tool is straightforward and one the field is already pursuing: training foundational UIPs on significantly more data may get us there even without further algorithmic or model improvements. If we manage to double or triple the discovery acceleration factor, the benefits of including ML in materials discovery workflows will significantly outweigh the setup costs.

We welcome further model submissions at
[{repo}]({repo}).

## Acknowledgments

Janosh Riebesell acknowledges support from the German Academic Scholarship Foundation ([Studienstiftung](https://wikipedia.org/wiki/Studienstiftung)).

A big thank you to

- Hai-Chen Wang and co-authors for creating and freely providing the WBM data set
- Jason Blake Gibson, Shyue Ping Ong, Chi Chen, Tian Xie, Bowen Deng, Peichen Zhong, Ekin Dogus Cubuk for helpful discussions
- Philipp Benner ([@pbenner](https://github.com/pbenner)) for [finding and reporting many bugs]({repo}/issues?q=is%3Aissue+author%3Apbenner+) in the data loading routines before the v1 release.

## Author Contributions

Janosh Riebesell: Methodology, Software, Data Curation, Training and Testing Models, Formal Analysis. Rhys Goodall: Conceptualization, Software, Formal Analysis. Anubhav Jain: Supervision. Kristin Persson: Supervision. Alpha Lee: Supervision.
