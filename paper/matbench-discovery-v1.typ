#import "template.typ": arkheion, arkheion-appendices, pdf-img, subfigure

#let citation = yaml(read("../citation.cff", encoding: none))
#let authors = (
  citation.authors.map(author => (
    name: author.given-names + " " + author.family-names,
    orcid: author.orcid,
    affiliation: author.affiliation,
    email: if author.keys().contains("email") { author.email } else { "" },
  ))
)

#show: arkheion.with(
  title: citation.title + ": A framework to evaluate machine learning crystal stability predictions",
  authors: authors,
  abstract: [
    The rapid adoption of machine learning in various scientific domains calls for the development of best practices and community agreed-upon benchmarking tasks and metrics. We present Matbench Discovery as an example evaluation framework for machine learning (ML) energy models, here applied as pre-filters to first-principles computed data in a high-throughput search for stable inorganic crystals.
    We address the disconnect between (i) thermodynamic stability and formation energy and (ii) retrospective and prospective benchmarking for materials discovery.
    Alongside this paper, we publish a Python package to aid with future model submissions and a growing online leaderboard with adaptive user-defined weighting of various performance metrics allowing researchers to prioritize the metrics they value most.
    To answer the question of which ML methodology performs best at materials discovery, our initial release includes random forests, graph neural networks (GNN), one-shot predictors, iterative Bayesian optimizers and foundational force fields (FFF).
    We highlight a misalignment between commonly used regression metrics and more task-relevant classification metrics for materials discovery.
    Accurate regressors are susceptible to unexpectedly high false-positive rates if those accurate predictions lie close to the decision boundary at 0 eV/atom above the convex hull.
    The benchmark results demonstrate that FFFs have advanced sufficiently to effectively and cheaply pre-screen thermodynamic stable hypothetical materials in future expansions of high-throughput materials databases.
  ],
  keywords: (
    "Materials Discovery",
    "Machine Learning",
    "Interatomic Potentials",
    "Foundation Models",
    "Density Functional Theory",
    "Benchmarking",
    "Convex Hull of Stability",
  ),
  date: "Feb 27, 2025", // date of acceptance at Nature Machine Intelligence
)

= Introduction
<sec:introduction>
The challenge of evaluating, benchmarking and then applying the rapid evolution of machine learning models is common across scientific domains. Specifically, the lack of agreed-upon tasks and data sets can obscure the performance of the model, making comparisons difficult. Material science is one such domain, where in the last decade, the number of ML publications and associated models have increased dramatically. Similar to other domains, such as drug discovery and protein design, the ultimate success is often associated with the discovery of a new material with specific functionality. In the combinatorial sense, material science can be viewed as an optimization problem of mixing and arranging different atoms with a merit function that captures the complex range of properties that emerge. To date, $~ 10^5$ combinations have been tested experimentally @bergerhoff_inorganic_1983@belsky_new_2002 , $~ 10^7$ have been simulated@jain_commentary_2013@saal_materials_2013@curtarolo_aflow_2012@draxl_nomad_2018@schmidt_improving_2024, and upwards of $~ 10^10$ possible quaternary materials are allowed by electronegativity and charge-balancing rules @davies_computational_2016. The space of quinternaries and higher is even less explored, leaving vast numbers of potentially useful materials to be discovered. The discovery of new materials is a key driver of technological progress and lies on the path to more efficient solar cells, lighter and longer-lived batteries, smaller and more efficient transistor gates just to name a few. In light of our sustainability goals, these advances cannot come fast enough. Any speed-up new discovery methods might yield should be leveraged to their fullest extent.

Computational materials discovery continues to present significant challenges despite advances in theory and methodology. The process typically requires performing extensive high-throughput calculations, which can be computationally intensive and time-consuming. Moreover, the complex relationship between structure and properties means that finding materials with desired characteristics often remains more art than science. Machine learning approaches offer promising alternatives by efficiently identifying patterns within large datasets. These methods excel at handling multidimensional data, balancing multiple optimization objectives @riebesell_pushing_2024, quantifying prediction uncertainty @borg_quantifying_2022@goodall_rapid_2022@zhu_fast_2023, and extracting meaningful information from sparse or noisy data @depeweg_decomposition_2017@bartel_critical_2020. These capabilities make ML particularly valuable as a complementary tool to traditional computational methods in materials science.

In particular, we focus on the role of ML to accelerate the use of Kohn-Sham density-functional theory (DFT) in the materials discovery pipeline. In comparison to other simulation frameworks DFT offers a compelling compromise between fidelity and cost that has seen it adopted as a workhorse method by the computational material science community. The great strengths of DFT as a methodology have led it to demand up to 45% of core hours at the UK-based Archer2 Tier 1 supercomputer @montanari_goldilocks_2024 and over 70% allocation time in materials science sector at NERSC @griffin_computational_2022@austin_nerscworkload_2018. This heavy resource requirement drives demand for ways to reduce or alleviate its computational burden, such as efficiency improvements or substitution from ML approaches.

While typically exhibiting lower accuracy and reliability, ML models produce results significantly faster---by orders of magnitude---than ab-initio simulations. This speed advantage positions them ideally for high-throughput (HT) screening campaigns, where they can act as efficient pre-filters for computationally demanding, higher-fidelity methods like DFT. The pioneering work of Behler and Parrinello @behler_generalized_2007 demonstrated the use of neural networks to learn the DFT potential energy surface (PES). This breakthrough spurred rapid advancements and extensive efforts to train increasingly sophisticated ML models on available PES data. Early applications often involved deploying these models as interatomic potentials (or force fields) focused on specific materials, a process necessitating the creation of bespoke training datasets for each system under investigation @bartok_machine_2018@deringer_generalpurpose_2020. As larger and more diverse datasets have emerged from initiatives like the Materials Project (MP) @jain_commentary_2013, AFLOW @curtarolo_aflow_2012 or the Open Quantum Materials Database (OQMD) @saal_materials_2013, researchers have begun to train so-called universal models that cover 90 or more of the most-application relevant elements in the periodic table. This opens up the prospect of ML-guided materials discovery to increase the hit rate of stable crystals and speed up DFT- and expert-driven searches.

Progress in ML for materials is often measured according to performance on standard benchmark datasets. As ML models have grown in complexity and applicability, benchmark datasets need to grow with them to accurately measure their usefulness. However, due to the rapid pace of the field and the variety of possible approaches for framing the discovery problem, no large-scale benchmark yet exists for measuring the ability of ML to accelerate materials discovery. As a result, it is unclear which methodologies or models are best suited for this task. The materials community has explored several approaches for computational discovery, including coordinate-free predictors that operate without requiring precise atomic positions @goodall_rapid_2022, sequential optimization methods based on Bayesian principles @zuo_accelerating_2021, and physics-informed interatomic potentials with universal element coverage @chen_universal_2022@deng_chgnet_2023@batatia_mace_2023. While each approach has demonstrated success in specific contexts, systematic comparison across methodologies has been lacking, preventing clear identification of optimal approaches for materials discovery at scale. Our work aims to identify the state-of-the-art (SOTA) model by proposing an evaluation framework that closely simulates a real-world discovery campaign guided by ML models. Our analysis reveals that FFFs surpass all other methodologies we evaluated in terms of both accuracy and robustness.

We hope that creation of benchmarks following this framework create a pathway through which interdisciplinary researchers with limited material science backgrounds can contribute usefully to model architecture and methodology development on a relevant task and thereby aid progress in Material Science. This work expands on initial research conducted in the first author's PhD thesis @riebesell_machine_2024.

= Evaluation Framework for Materials Discovery
<sec:evaluation-framework-for-materials-discovery>
This work proposes a benchmark task designed to address four fundamental challenges that we believe are essential to justify the effort of experimentally validating machine learning predictions:

+ #strong[Prospective Benchmarking]: Idealized and overly simplified benchmarks may not adequately capture the challenges encountered in real-world applications. This disconnect can arise from selecting inappropriate targets @bartel_critical_2020 or using unrepresentative data splits @wu_moleculenet_2018@kpanou_robustness_2021. For small datasets of materials properties, 'Leave-Out' data splitting strategies are often used to assess model performance @meredig_can_2018@cubuk_screening_2019@zahrt_cautionary_2020. However, in our target domain large quantities of diverse data ($~ 10^5$) are available and hence retrospective splitting strategies predicated on clustering can end up testing artificial or unrepresentative use cases. This encourages using new sources of prospectively generated test data to understand application performance. Adopting this principle, the intended discovery workflow should be used to generate the test data leading to a substantial but realistic covariate shift between the training and test distributions that gives a much better indicator likely of performance on additional application of the same discovery workflow.

+ #strong[Relevant Targets:] For materials discovery, high-throughput DFT formation energies are widely used as regression targets but don't directly indicate thermodynamic stability or synthesizability. The true stability of a material depends on its energetic competition with other phases in the same chemical system, quantified by the distance to the convex hull phase diagram. This distance serves as the primary indicator of (meta-)stability under standard conditions @sun_thermodynamic_2016, making it a more suitable target despite other factors like kinetic and entropic stabilization that influence real-world stability but are more expensive to simulate, especially at scale.

  Additionally, ML models that require relaxed structures as input create a circular dependency with the DFT calculations they're meant to accelerate, reducing their practical utility for discovery.

+ #strong[Informative Metrics]: Global metrics such as $upright("MAE")$, $upright("RMSE")$ and $R^2$ may provide practitioners with misleading confidence regarding model reliability. Even models with strong regression performance can produce unexpectedly high rates of false-positive predictions when nominally accurate estimates fall near decision boundaries, resulting in substantial opportunity costs through wasted laboratory resources and time. Consequently, models should be evaluated based on their ability to facilitate correct decision-making patterns rather than regression accuracy alone. One effective approach is to define selection criteria and assess regression models primarily by their classification performance.

+ #strong[Scalability]: Future materials discovery efforts are likely to target broad chemical spaces and large data regimes. Small benchmarks can lack chemical diversity, and obfuscate poor scaling relations or weak out-of-distribution (OOD) performance. For instance, random forests achieve excellent performance on small datasets but are typically outperformed by neural networks on large datasets due to the benefits of representation learning @goodall_predicting_2020. Whilst we propose that large training sets are necessary to adequately differentiate the ability of models to learn in the larger data regime, given the enormous size of the configurational space of materials yet to be explored we also propose that it is important to construct a task where the test set is larger than the training set to mimic true deployment at scale. No other inorganic materials benchmarks test the prospects of large-scale deployment in this manner.

We highlight two specific benchmarking efforts that have partially addressed the above challenges: Matbench @dunn_benchmarking_2020 and the Open Catalyst Project (OCP) @chanussot_open_2021. Other valuable efforts such as MatSciML @lee_foundation_2023 and JARVIS-Leaderboard @choudhary_jarvisleaderboard_2024 aggregate a wide variety of material science related benchmark tasks, including from Matbench and OCP, but do not introduce distinct benchmarking design patterns to those seen in Matbench or the OCP.

By providing a standardized collection of 13 datasets ranging in size from $~$300 to $~$132,000 samples from both DFT and experimental sources, Matbench addresses the scalability challenge, highlighting how model performance changes as a function of data regime. Matbench helped focus the field of ML for materials, increase comparability across papers and provide a quantitative measure of progress in the field. Importantly, all tasks were exclusively concerned with the properties of known materials. We believe a task that simulates a materials discovery campaign by requiring materials stability prediction from unrelaxed structures to be a missing piece here.

OCP is a large-scale initiative aimed at discovering substrate-adsorbate combinations that can catalyze critical industrial reactions, transforming these adsorbates into more useful products. The OCP has released two datasets thus far, OCP20 @chanussot_open_2021 and OCP22 @tran_open_2022 for training and benchmarking ML models. OCP certainly addressed challenge 1 of closely mimicking a real-world problem. They have recently shown that despite not reaching their target accuracy to entirely replace DFT, using ML in conjunction with confirmatory DFT calculations dramatically speeds up their combinatorial screening workflow @lan_adsorbml_2023. The team behind the OCP has a second initiative targeting materials for direct air capture called OpenDAC that has shared the ODAC23 dataset @sriram_open_2023. The OpenDAC benchmark is setup identically to the OCP.

We believe that addressing these four challenges will result in benchmarks that enable future ML-guided discovery efforts to confidently select appropriate models and methodologies for the expansion of computational materials databases.

#figure(
  pdf-img("figs/overview.pdf"),
  caption: [
    An overview of how data is used in Matbench-Discovery. a) shows a conventional prototype-based discovery workflow where different elemental assignments to the sites in a known prototype are used to create a candidate structure. This candidate is relaxed using DFT to arrive at a relaxed structure that can be compared against a reference convex hull. This sort of workflow was used to construct the WBM data set. b) highlights that databases such as the Materials Project provide a rich set of data which different academic groups have used to explore different types of models. While prior work tended to focus on individual modalities, our framework enables consistent model comparisons across modalities. c) shows the proposed test evaluation framework where the end user takes a machine learning model and uses it to predict a relaxed energy given an initial structure (IS2RE). This energy is then used to make a prediction as to whether the material will be stable or unstable with respect to a reference convex hull. From an applications perspective, this classification performance is better aligned with intended use cases in screening workflows.
  ],
)<fig:overview>

= Results
<sec:results>
@tab:metrics-table-uniq-protos shows performance metrics for all models included in the initial release of Matbench Discovery reported on the unique protostructure subset. EquiformerV2 + DeNS achieved the highest performance in ML-guided materials discovery, surpassing all other models across the 9 reported metrics. When computing metrics in the presence of missing values or obviously pathological predictions (error of 5 eV/atom or greater) we assign the dummy regression values and a negative classification prediction to these points. The discovery acceleration factor (DAF) quantifies how many times more effective a model is at finding stable structures compared to random selection from the test set. Formally the DAF is the ratio of the precision to the prevalence. The maximum possible DAF is the inverse of the prevalence which on our dataset is $\( upright("33k") \/ upright("215k") \)^(- 1) approx 6.5$. Thus the current SOTA of 5.04 achieved by EquiformerV2 + DeNS leaves room for improvement. However, evaluating each model on the subset of 10k materials each model ranks as being most stable (see @tab:metrics-table-first-10k, we see an impressive DAF of 6.33 for EquiformerV2 + DeNS which is approaching optimal performance for this task.

A significant performance gap emerges between models predicting energy directly from unrelaxed inputs (MEGNet, Wrenformer, CGCNN, CGCNN+P, ALIGNN, Voronoi RF) and FFFs, which leverage force and stress data to emulate DFT relaxation for final energy prediction. While the energy-only models exhibit surprisingly strong classification metrics (F1, DAF), their regression performance ($R^2$, RMSE) is considerably poorer. Notably, only ALIGNN, BOWSR, and CGCNN+P among the energy-only models achieve a positive coefficient of determination ($R^2$). Negative $R^2$ means model predictions explain the observed variation in the data less than simply predicting the test set mean. In other words, these models are not predictive in a global sense (across the full dataset range). Nevertheless, models with negative $R^2$ may still show predictive capability for materials far from the stability threshold (i.e., in the distribution tails). Their performance suffers most near the 0~eV/atom stability threshold, the region with the highest concentration of materials. This illustrates a limitation of using $R^2$ alone to evaluate models for classification tasks such as stability prediction.

#figure(
  pdf-img("figs/cumulative-precision-recall-only-compliant.pdf"),
  caption: [
    Precision and recall as a function of the number of model predictions validated. A typical discovery campaign will rank hypothetical materials by model-predicted hull distance from most to least stable and validate the most stable predictions first. A higher fraction of correct stable predictions corresponds to higher precision and fewer stable materials overlooked correspond to higher recall. Precision is calculated based only on the selected materials up to that point, whilst the cumulative recall depends on knowing the total number of positives upfront. Models like eqV2 S DeNS and Orb MPtrj perform better for exhaustive discovery campaigns (screening a higher share of the candidate pool), others like CHGNet do better when validating a smaller percentage of the materials predicted to be most stable. FFFs offer significantly improved precision on shorter campaigns of $~$20k or less materials validated, as they are less prone to false positive predictions among highly stable materials.
  ],
)<fig:cumulative-precision-recall>

The reason CGCNN+P achieves better regression metrics than CGCNN but is still worse as a classifier becomes apparent from @fig:hist-clf-pred-hull-dist-models by noting that the CGCNN+P histogram is more sharply peaked at the 0 hull distance stability threshold. This causes even small errors in the predicted convex hull distance to be large enough to invert a classification. Again, this is evidence to choose carefully which metrics to optimize. Regression metrics are far more prevalent when evaluating energy predictions. However, our benchmark treats energy predictions as merely means to an end to classify compound stability. Improvements in regression accuracy are of limited use to materials discovery in their own right unless they also improve classification accuracy. Our results demonstrate that this is not a given.

@fig:cumulative-precision-recall has models rank materials by model-predicted hull distance from most to least stable; materials furthest below the known hull at the top, materials right on the hull at the bottom. For each model, we iterate through that list and calculate at each step the precision and recall of correctly identified stable materials. This simulates exactly how these models would be used in a prospective materials discovery campaign and reveals how a model's performance changes as a function of the discovery campaign length. As a practitioner, you typically have a certain amount of resources available to validate model predictions. These curves allow you to read off the best model given these constraints. For instance, plotting the results in this manner shows that CHGNet initially achieves higher precision, than models such as EquiformerV2+DeNS, ORB MPTrj, SevenNet, and MACE that report higher precision across the whole test set.

In @fig:cumulative-precision-recall each line terminates when the model believes there are no more materials in the WBM test set below the MP convex hull. The dashed vertical line shows the actual number of stable structures in our test set. All models are biased towards stability to some degree as they all overestimate this number, most of all BOWSR by 133%. This overestimation primarily affects exhaustive discovery campaigns aiming to validate all materials predicted as stable. In practice, campaigns are often resource-limited (e.g., to 10,000 DFT relaxations). By ranking candidates by predicted stability and validating only the top fraction dictated by the budget, the higher concentration of false positives typically found among less stable predictions is avoided without diminishing the campaign's effective discovery rate (see @tab:metrics-table-first-10k where even the DAF of the worst performing model benchmarked, Voronoi RF, jumps from 1.58 to 2.49).

The diagonal 'Optimal Recall' line on the recall plot in @fig:cumulative-precision-recall would be achieved if a model never made a false negative prediction and stopped predicting stable crystals exactly when the true number of stable materials is reached. Examining the FFF models, we find that they all achieve similar recall values, ranging from approximately 0.75 to 0.86.. This is substantially smaller than the variation we see in the precision for the same models $~$0.44-0.77. Inspecting the overlap, we find that the intersection of the models' correct agreements accounts for a precision of only 0.57 within the $~$0.75-0.86 range, with just 0.04 of the examples where all models are wrong simultaneously. These results indicate that the models are making meaningfully different predictions.

Examining the precision plot in @fig:cumulative-precision-recall, we observe that the energy-only models exhibit a much more pronounced drop in their precision early-on, falling to 0.6 or less in the first 5k screened materials. Many of these models (all except BOWSR, Wrenformer and Voronoi RF) display an interesting hook shape in their cumulative precision, recovering again slightly in the middle of the simulated campaign between ~5k and up to ~30k before dropping again until the end.

#figure(
  pdf-img("figs/rolling-mae-vs-hull-dist-models-only-compliant.pdf"),
  caption: [
    Universal potentials are more reliable classifiers because they exit the red triangle earliest. The lines represent rolling mean absolute error (MAE) on the WBM test set as a function of distance to the MP training set convex hull. The red 'triangle of peril' indicates regions where the mean error exceeds the distance to the stability threshold (0 eV). Within this triangle, models are more likely to misclassify materials as the errors can flip classifications. Earlier exit from the triangle correlates with fewer false positives (right side) or false negatives (left side). The width of the 'rolling window' indicates the range over which prediction errors were averaged.
  ],
)<fig:rolling-mae-vs-hull-dist-models>

@fig:rolling-mae-vs-hull-dist-models provides a visual representation of the reliability of different models as a function of a material's DFT distance to the MP convex hull. The lines show the rolling mean absolute error (RMAE) of model-predicted hull distances vs DFT. The red-shaded area, which we coin the 'triangle of peril', emphasizes the zone where the average model error surpasses the distance to the stability threshold at 0 eV. As long as the rolling MAE remains within this triangle, the model is highly susceptible to misclassifying structures. The average error in this region is larger than the distance to the classification threshold at 0, consequently in cases where the error points towards the stability threshold it would be large enough to flip a correct classification into an incorrect one. Inside this region, the average error magnitude surpasses the distance to the classification threshold at 0 eV. Consequently, when errors point toward the stability boundary, they are sufficiently large to potentially reverse a correct classification. The faster a model's error curve exits the triangle on the left side (representing negative DFT hull distances), the lower its tendency to mistakenly classify stable structures as unstable, thereby reducing false negatives. Exiting promptly on the right side (positive DFT hull distances) correlates with a decreased probability of predicting unstable structures as stable, resulting in fewer false positives.

Models generally exhibit lower rolling errors towards the left edge of the plot compared to the right edge. This imbalance indicates a greater inclination towards false positive predictions than false negative ones. Put differently, all models are less prone to predicting a material at −0.2~eV/atom DFT hull distance as unstable than they are to predicting a material at 0.2~eV/atom DFT hull distance as stable. From a practical perspective, this is undesirable because the opportunity cost associated with validating an incorrectly predicted stable material (a false positive) is typically much higher than that of missing a genuinely stable one (a false negative). We hypothesize that this error asymmetry arises from the MP training set's uncharacteristically high proportion of stable materials, causing statistical models trained on it to be biased towards assigning low energies even to high-energy atomic arrangements. Training on datasets with more high-energy structures, such as Alexandria @schmidt_improving_2024 and OMat24 @barroso-luque_open_2024, would be expected to improve performance by balancing out this source of bias.

= Discussion
<sec:discussion>
We have demonstrated the effectiveness of ML-based triage in high-throughput materials discovery and posit that the benefits of including ML in discovery workflows now clearly outweigh the costs. @tab:metrics-table-uniq-protos shows that in a realistic benchmark scenario that several models achieve a discovery acceleration greater than 2.5 across the whole dataset and up to 6 when considering only the 10k most stable predictions from each model (@tab:metrics-table-first-10k). Initially, the most promising ML methodology for accelerating high-throughput discovery was uncertain. Our results reveal a distinct advantage for FFFs regarding both accuracy and extrapolation performance. Incorporating force information allows FFFs to better simulate the relaxation pathway towards the DFT-relaxed structure, enabling a more accurate final energy determination.

Ranked best-to-worst by their test set F1 score on thermodynamic stability prediction, we find EquiformerV2 + DeNS \> Orb \> SevenNet \> MACE \> CHGNet \> M3GNet \> ALIGNN \> MEGNet \> CGCNN \> CGCNN+P \> Wrenformer \> BOWSR \> Voronoi fingerprint random forest. The top models are FFFs which we establish to be the best methodology for ML-guided materials discovery, achieving F1 scores of 0.57-0.82 for crystal stability classification and discovery acceleration factors (DAF) of up to 6x on the first 10k most stable predictions compared to dummy selection.

As the convex hull becomes more comprehensively sampled through future discoveries, the fraction of unknown stable structures will naturally decline. This will lead to less enriched test sets and, consequently, more challenging and discriminative discovery benchmarks. However, the discovery task framed here addresses only a limited subset of potential FFF applications. We believe that additional benchmarks are essential to effectively guide FFF development. These efforts should prioritize task-based evaluation frameworks that address the four critical challenges we identify for narrowing the deployment gap: adopting prospective rather than retrospective benchmarking, tackling relevant targets, using informative metrics, and scalability.

Looking ahead, the consistently linear log-log learning curves observed in related literature @vonlilienfeld_retrospective_2020 suggest that further decreases in the error of FFFs can be readily unlocked with increased training data. This has be borne out in the scaling results of GNoME @merchant_scaling_2023, MatterSim @yang_mattersim_2024, Alexandria @schmidt_improving_2024, and OMat24 @barroso-luque_open_2024 which all show improvements in performance when training on much larger datasets. To realize the full potential of scaling these models, future efforts should deploy their resources to generate large quantities of higher-than-PBE fidelity training data. The quality of a FFF model is circumscribed by the quality and level of theory of its training data.

Beyond simply predicting thermodynamic stability at zero Kelvin, future models will need to understand and predict material properties under varying environmental conditions, such as finite temperature and pressure, to aid in materials discovery. In this context, temperature-dependent dynamical properties constitute an area ripe for interatomic potentials. Another key open question is how effectively these models can contribute to the computational prediction of synthesis pathways. Many current methods for predicting reaction pathways employ heuristic rules to manage the considerable complexity introduced by metastability, in addition to relying on conventional ground-state ab initio data @mcdermott_graphbased_2021@aykol_rational_2021@wen_chemical_2023. These algorithms will massively benefit from more efficient estimates of reaction energy barriers @yuan_analytical_2024 and non-crystalline, out-of-equilibrium materials @aykol_thermodynamic_2018, opening up a whole new field to machine learning accelerated inquiry.

= Methods
<sec:methods>
== Matbench Discovery Framework
<sec:matbench-discovery>
As first presented in the first author's PhD thesis @riebesell_machine_2024, we propose an evaluation framework that places no constraints on the type of data a model is trained on as long as it would be available to a practitioner conducting a real materials discovery campaign. This means that for the high-throughput DFT data considered, any subset of the energies, forces, stresses or any other properties that can be routinely extracted from DFT calculations, such as magnetic moments, are valid training targets. All of these would be available to a practitioner performing a real materials discovery campaign and hence are permitted for training any model submission. We only enforce that at test time, all models must make predictions on the convex hull distance of the relaxed structure with only the unrelaxed structure as input. This setup avoids circularity in the discovery process, as unrelaxed structures can be cheaply enumerated through elemental substitution methodologies and do not contain information inaccessible in a prospective discovery campaign. @fig:overview provides a visual overview of design choices.

The convex hull distance of a relaxed structure is chosen as the measure of its thermodynamic stability, rather than the formation energy, as it informs the decision on whether to pursue a potential candidate crystal. This decision was also motivated by @bartel_critical_2020 who found that even composition-only models are capable of predicting DFT formation energies with useful accuracy. However, when tasking those same models with predicting decomposition enthalpy, performance deteriorated sharply. This insight exposes how ML models are much less useful than DFT for discovering new inorganic solids than would be expected given their low prediction errors for formation energies due to the impact of random as opposed to systematic errors.

Standard practice in ML benchmarks is to hold all variables fixed -- most importantly the training data -- and vary only the model architecture to isolate architectural effects on the performance. We deliberately deviate from this practice due to diverging objectives from common ML benchmarks. Our goal is to identify the best methodology for accelerating materials discovery. What kind of training data a model can ingest is part of its methodology. Unlike energy-only models, FFFs benefit from the additional training data provided by the forces and stresses recorded in DFT relaxations. This allows them to learn a fundamentally higher fidelity model of the physical interactions between ions. That is a genuine advantage of the architecture and something any benchmark aiming to identify the optimal methodology for materials discovery must reflect. In light of this utilitarian perspective, our benchmark contains models trained on varying datasets, and any model that can intake more physical modalities from DFT calculations is a valid model for materials discovery.

We define the Materials Project (MP) @jain_commentary_2013 #link("https://docs.materialsproject.org/changes/database-versions")[v2022.10.28] database release as the maximum allowed training set for any compliant model submission. Models may train on the complete set of relaxation frames, or any subset thereof such as the final relaxed structures. Any subset of the energies, forces, and stresses are valid training targets. In addition any auxiliary tasks such as predicting electron densities, magnetic moments, site-partitioned charges, etc. that can be extracted from the output of the DFT calculations are allowed for multi-task learning @shoghi_molecules_2023. Our test set consists of the unrelaxed structures in the WBM dataset @wang_predicting_2021. Their target values are the PBE formation energies of the corresponding DFT-relaxed structures.

=== Limitations of this Framework
<sec:limitations>
Whilst the framework proposed here directly mimics a common computational materials discovery workflow, it is worth highlighting that there still exist significant limitations to these traditional computational workflows that can prevent the material candidates suggested by such a workflow from being able to be synthesized in practice. For example, high-throughput DFT calculations often use small unit cells which can lead to artificial orderings of atoms. The corresponding real material may be disordered due to entropic effects that cannot be captured in the zero-kelvin thermodynamic convex hull approximated by DFT @cheetham_artificial_2024.

Another issue is that, when considering small unit cells, the DFT relaxations may get trapped at dynamically unstable saddle points in the true PES. This failure can be detected by calculating the phonon spectra for materials predicted to be stable. However, the cost of doing so with DFT is often deemed prohibitive for high-throughput screening. The lack of information about the dynamic stability of nominally stable materials in the WBM test set prevents this work from considering this important criteria as an additional screening filter. However, recent progress in the development of FFFs suggests that ML approaches will soon provide sufficiently cheap approximations of these terms for high-throughput searches @batatia_foundation_2023@deng_overcoming_2024. As the task presented here begins to saturate, we believe that future discovery benchmarks should extend upon the framework proposed here to also incorporate criteria based on dynamic stability.

When training FFF models there is a competition between how well given models can fit the energies, forces, and stresses simultaneously. The metrics in the Matbench Discovery leaderboard are skewed towards energies and consequently FFF models trained with higher weighting on energies can achieve better metrics. We caution that optimizing hyperparameters purely to improve performance on this benchmark may have unintended consequences for models intended for general purpose use. Practitioners should also consider other involved evaluation frameworks that explore orthogonal use cases when developing model architectures. We highlight work from #cite(<pota_thermal_2024>, form: "prose") on thermal conductivity benchmarking, #cite(<fu_forces_2023>, form: "prose") on MD stability for molecular simulation, and #cite(<chiang_mlip_2025>, form: "prose") on modeling reactivity (hydrogen combustion) and asymptotic correctness (homonuclear diatomic energy curves) as complementary evaluation tasks for assessing the performance of FFF models.

We design the benchmark considering a positive label for classification as being on or below the convex hull of the MP training set. An alternative formulation would be to say that materials in WBM that are below the MP convex hull but do not sit on the combined MP+WBM convex hull are negatives. The issue with such a design is that it involves unstable evaluation metrics. If we consider the performance against the final combined convex hull rather than the initial MP convex hull, then each additional sample considered can retroactively change whether or not a previous candidate would be labeled as a success as it may no-longer sit on the hull. Since constructing the convex hull is computationally expensive, this path dependence makes it impractical to evaluate cumulative precision metrics (see @fig:cumulative-precision-recall). The chosen setup does increase the number of positive labels and could consequently be interpreted as overestimating the performance. This overestimation decreases as the convex hull becomes better sampled. Future benchmarks building on this work could make use of combination of MP+WBM to control this artifact. An alternative framework could report metrics for each WBM batch in turn and retrain between batches, this approach was undesirable here as it increases the cost of submission five-fold and introduces many complexities, for example, should each model only retrain on candidates it believed to be positive, that would make fair comparison harder.

== Datasets
<sec:datasets>
=== Materials Project Training Set
<sec:materials-project-training-set>
The Materials Project is a widely-used database of inorganic materials properties that have been calculated using high-throughput ab-initio methods. At the time of writing, the Materials Project database @jain_commentary_2013 has grown to $~$154 k crystals, covering diverse chemistries and providing relaxed and initial structures as well as the relaxation trajectory for every entry.

Our benchmark defines the training set as all data available from the #link("https://docs.materialsproject.org/changes/database-versions#v2022.10.28")[v2022.10.28 MP release]. We recorded a snapshot of energies, forces, stresses and magnetic moments for all MP ionic steps on 2023-03-15 as the canonical training set for Matbench Discovery, and provide convenience functions through our #link("https://pypi.org/project/matbench-discovery")[Python package] for easily feeding that data into future model submissions to our benchmark.

Flexibility in specifying the dataset allows authors to experiment with and fully exploit the available data. This choice is motivated by two factors. First, it may seem that models trained on multiple snapshots containing energies, forces, and stresses receive more training data than models trained only on the energies of relaxed structures. However, the critical factor is that all this additional data was generated as a byproduct of the workflow to produce relaxed structures. Consequently, all models are being trained using data acquired at the same overall cost. If some architectures or approaches can leverage more of this byproduct data to make improved predictions this is a fair comparison between the two models. This approach diverges philosophically from other benchmarks such as the OCP and Matbench where it has been more common to subcategorize different models and look at multiple alternative tasks (e.g. composition-only vs structure-available in Matbench or IS2RS, IS2RE, S2EF in OCP) and do not make direct comparisons of this manner. Second, recent work in the space from @li_critical_2023@li_exploiting_2023 has claimed that much of the data in large databases like MP are redundant and that models can be trained more effectively by taking a subset of these large data pools. From a systems-level perspective, identifying innovative cleaning or active-learning strategies to make better use of available data may be as crucial as architectural improvements, as both can similarly enhance performance, especially given the prevalence of errors in high-throughput DFT. Consequently, such strategies where they lead to improved performance should be able to be recognized within the benchmark. We encourage submissions to submit ablation studies showing how different system-level choices affect performance. Another example of a system-level choice that may impact performance is the choice of optimizer, for example FIRE @bitzek_structural_2006 vs L-BFGS, in the relaxation when using FFF models.

We highlight several example datasets that are valid within the rules of the benchmark that take advantage of these freedoms. The first is the MP-crystals-2019.4.1 dataset @chen_graph_2019 which is a subset of 133,420 crystals and their formation energies that form a subset of the v2021.02.08 MP release. The MP-crystals-2022.10.28 dataset is introduced with this work comprising a set of 154,719 structures and their formation energies drawn from the v2021.02.08 MP release. The next is the MPF.2021.2.8 dataset @chen_universal_2022 curated to train the M3GNet model which takes a subset of 62,783 materials from the v2021.02.08 MP release. The curators of the MPF.2021.2.8 dataset down-sampled the v2021.02.08 release significantly to select a subset of calculations that they believed to be most self-consistent. Rather than taking every ionic step from the relaxation trajectory this dataset opts to select only the initial, final and one intermediate structure for each material to avoid biasing the dataset towards examples where more ionic steps were needed to relax the structure. Consequently the dataset consists of 188,349 structures. The MPF.2021.2.8 is a proper subset of the training data as no materials were deprecated between the v2021.02.08 and v2022.10.28 database releases. The final dataset we highlight, with which several of the FFF models have been trained, is the MPtrj dataset @deng_chgnet_2023. This dataset was curated from the earlier v2021.11.10 MP release. The MPtrj dataset is a proper subset of the allowed training data but several potentially anomalous examples from within MP were cleaned out of the dataset before the frames were subsampled to remove redundant frames. It is worth noting that the v2022.10.28 release contains a small number of additional Perovskite structures not found in MPtrj that could be added to the training set within the scope of the benchmark.

We note that the v2023.11.1 deprecated a large number of calculations so data queried from subsequent database releases is not considered valid for this benchmark.

=== WBM Test Set
<sec:wbm-test-set>
The WBM dataset @wang_predicting_2021 consists of 257,487 structures generated via chemical similarity-based elemental substitution of MP source structures followed by DFT relaxation and calculating each crystal's convex hull distance. The element substitutions applied to a given source structure were determined by random sampling according to the weights in a chemical similarity matrix data-mined from the ICSD @glawe_optimal_2016.

The WBM authors performed 5 iterations of this substitution process (we refer to these steps as batches). After each step, the newly generated structures found to be thermodynamically stable after DFT relaxation flow back into the source pool to partake in the next round of substitution. This split of the data into batches of increasing substitution count is a unique and compelling feature of the test set as it allows out-of-distribution (OOD) testing by examining whether model performance degrades for later batches. A higher number of elemental substitutions on average carries the structure further away from the region of material space covered by the MP training set (see @fig:rolling-mae-vs-hull-dist-wbm-batches-models for details). Whilst this batch information makes the WBM dataset an exceptionally useful data source for examining the extrapolation performance of ML models, we look primarily at metrics that consider all batches as a single test set.

In order to control for the potential adverse effects of leakage between the MP training set and the WBM test set, we cleaned the WBM test set based on protostructure matching. We refer to the combination of a materials prototype and the elemental assignment of its wyckoff positions as a protostructure following @parackal_identifying_2024. First we removed 524 pathological structures in WBM based on formation energies being larger than 5 eV/atom or smaller than -5 eV/atom. We then removed from the WBM test set all examples where the final protostructure of a WBM material matched the final protostructure of an MP material. In total 11,175 materials were cleaned using this filter. We further removed all duplicated protostructures within WBM, keeping the lowest energy structure in each instance, leaving 215,488 structures in the unique prototype test set.

Throughout this work, we define stability as being on or below the convex hull of the MP training set ($E_(upright("MP hull dist")) lt.eq 0$). 32,942 out of 215,488 materials in WBM unique prototype test set satisfy this criterion. Of these, $~$33k are unique prototypes, meaning they have no matching structure prototype in MP nor another higher-energy duplicate prototype in WBM. Our code treats the stability threshold as a dynamic parameter, allowing for future model comparisons at different thresholds. For initial analysis in this direction, see @fig:roc-models in the SI.

As WBM explores regions of materials space not well sampled by MP, many of the discovered materials that lie below MP's convex hull are not stable relative to each other. Of the $~$33k that lie below the MP convex hull less than half, or around $~$20k, remain on the joint MP+WBM convex hull. This observation suggests that many WBM structures are repeated samples in the same new chemical spaces. It also highlights a critical aspect of this benchmark in that we knowingly operate on an incomplete convex hull. Only current knowledge of competing points on the PES is accessible to a real discovery campaign and our metrics are designed to reflect this.

== Models
<sec:models>
#figure(
  pdf-img("figs/metrics-table-uniq-protos-only-compliant.pdf"),
  caption: [
    Metrics on the unique prototype test set.
  ],
)<tab:metrics-table-uniq-protos>
To test a wide variety of methodologies proposed for learning the potential energy landscape, our initial benchmark release includes 13 models. Next to each model's name we give the training targets that were used: E - Energy, F - Forces, S - Stresses and M - Magnetic moments. The subscripts G and D refer to whether gradient-based or direct prediction methods were used to obtain force and stress predictions.

+ #strong[EquiformerV2 + DeNS] @liao_equiformerv2_2024@liao_generalizing_2024 ($"EFS"_"D"$) - EquiformerV2 builds on the first Equiformer model @liao_equiformer_2023 by replacing the $"SO"\(3\)$ convolutions with eSCN convolutions @passaro_reducing_2023 as well as a range of additional tweaks to make better use of the ability to scale to higher $L_(upright("max"))$ using eSCN convolutions. EquiformerV2 uses direct force prediction rather than taking the forces as the derivative of the energy predictions for computational efficiency. Here we take the pre-trained `eqV2 S DeNS` @barroso-luque_open_2024 trained on the MPTrj dataset. This model in addition to supervised training on energies, forces, and stresses makes use of a auxiliary training task based on de-noising non-equilibrium structures @liao_generalizing_2024. We refer to this model as `EquiformerV2 + DeNS` in the text and `eqV2 S DeNS` in plots.

+ #strong[Orb] @neumann_orb_2024 ($"EFS"_"D"$) - Orb is a lightweight model architecture developed to scale well for the simulation of large systems such as metal organic frameworks. Rather than constructing an architecture that is equivariant by default, Orb instead makes use of data augmentation during training to achieve approximate equivariance. This simplifies the architecture allowing for faster inference. We report results for the `orb-mptrj-only-v2` model which was pre-trained using a diffusion-based task on MPTrj before supervised training on the energies, forces and stresses in MPTrj. For simplicity we refer to this model as `ORB MPTrj`.

+ #strong[SevenNet] @park_scalable_2024 ($"EFS"_"G"$) - SevenNet emerged from an effort to improve the performance of message passing neural networks @gilmer_neural_2017 when used to conduct large scale simulations involving that benefit from parallelism via spatial decomposition. Here we use the pre-trained `SevenNet-0_11July2024` trained on the MPTrj dataset. The SevenNet-0 model is an equivariant architecture based on a NeqFFF @batzner_e3equivariant_2022 architecture that mostly adopts the GNoME @merchant_scaling_2023 hyper-parameters. SevenNet-0 differs from NeqFFF and GNoME by replacing the tensor product in the self-connection layer with a linear layer applied directly to the node features, this reduces the number of parameters from 16.24 million in GNoME to 0.84 million for SevenNet-0. For simplicity we refer to this model as `SevenNet`.

+ #strong[MACE] @batatia_mace_2023 ($"EFS"_"G"$) - MACE builds upon the recent advances @thomas_tensor_2018@batzner_e3equivariant_2022 in equivariant neural network architectures by proposing an approach to computing high-body-order features efficiently via Atomic Cluster Expansion @drautz_atomic_2019. Unlike the other FFF models considered MACE was primarily developed for molecular dynamics of single material systems and not the universal use case studied here. The authors trained MACE on the MPTrj dataset, these models have been shared under the name `MACE-MP-0` @batatia_foundation_2023 and we report results for the `2023-12-03` version commonly called `MACE-MP-0 (medium)`. For simplicity we refer to this model as `MACE`.

+ #strong[CHGNet] @deng_chgnet_2023 ($"EFS"_"G"$M) - CHGNet is a FFF for charge-informed atomistic modeling. Its distinguishing feature is that it was trained to predict magnetic moments on top of energies, forces and stresses in the MPtrj dataset (which was prepared for the purposes of training CHGNet). By modeling magnetic moments, CHGNet learns to accurately represent the orbital occupancy of electrons which allows it to predict both atomic and electronic degrees of freedom. We make use of the pre-trained `v0.3.0` CHGNet model from @deng_chgnet_2023.

+ #strong[M3GNet] @chen_universal_2022 ($"EFS"_"G"$) - M3GNet is a GNN-based FFF for materials trained on up to 3-body interactions in the initial, middle and final frame of MP DFT relaxations. The model takes the unrelaxed input and emulates structure relaxation before predicting energy for the pseudo-relaxed structure. We make use of the pre-trained `v2022.9.20` M3GNet model from @chen_universal_2022 trained on the compliant MPF.2021.2.8 dataset.

+ #strong[ALIGNN] @choudhary_atomistic_2021 (E) - The Atomistic Line Graph Neural Network (ALIGNN) is a message passing GNN architecture that takes as input both the interatomic bond graph and a line graph corresponding to 3-body bond angles. The ALIGNN architecture involves a global pooling operation which means that it is ill-suited to force-field applications. To address this the ALIGNN-FF model was later introduced without global pooling @choudhary_unified_2023. We trained ALIGNN on the MP-crystals-2022.10.28 dataset for this benchmark.

+ #strong[MEGNet] @chen_graph_2019 (E) - MatErials Graph Network is another GNN-based architecture that also updates a set of edge and global features (like pressure and temperature) in its message passing operation. This work showed that learned element embeddings encode periodic chemical trends and can be transfer-learned from large datasets (formation energies) to predictions on small data properties (band gaps, elastic moduli). We make use of the pre-trained `Eform_MP_2019` MEGNet model trained on the compliant MP-crystals-2019.4.1 dataset.

+ #strong[CGCNN] @xie_crystal_2018 (E) - The Crystal Graph Convolutional Neural Network (CGCNN) was the first neural network model to directly learn 8 different DFT-computed material properties from a graph representing the atoms and bonds in a periodic crystal. CGCNN was among the first to show that just like in other areas of ML, given large enough training sets, neural networks can learn embeddings that outperform human-engineered structure features directly from the data. We trained an ensemble of 10 CGCNN models on the MP-crystals-2022.10.28 dataset for this benchmark.

+ #strong[CGCNN+P] @gibson_dataaugmentation_2022 (E) - This work proposes simple, physically motivated structure perturbations to augment stock CGCNN's training data of relaxed structures with structures resembling unrelaxed ones but mapped to the same DFT final energy. Here we chose $P = 5$, meaning the training set is augmented with 5 random perturbations of each relaxed MP structure mapped to the same target energy. In contrast to all other structure-based GNNs considered in this benchmark, CGCNN+P is not attempting to learn the Born-Oppenheimer potential energy surface. The model is instead taught the PES as a step-function that maps each valley to its local minimum. The idea is that during testing on unrelaxed structures, the model will predict the energy of the nearest basin in the PES. The authors confirm this by demonstrating a lowering of the energy error on unrelaxed structures. We trained an ensemble of 10 CGCNN+P models on the MP-crystals-2022.10.28 dataset for this benchmark.

+ #strong[Wrenformer] (E) - For this benchmark, we introduce Wrenformer which is a variation on the coordinate-free Wren model @goodall_rapid_2022 constructed using standard QKV-self-attention blocks @vaswani_attention_2017 in place of message-passing layers. This architectural adaptation reduces the memory usage allowing the architecture to scale to structures with greater than 16 Wyckoff positions. Like its predecessor, Wrenformer is a fast coordinate-free model aimed at accelerating screening campaigns where even the unrelaxed structure is a priori unknown @parackal_identifying_2024. The key idea is that by training on the coordinate anonymized Wyckoff positions (symmetry-related positions in the crystal structure), the model learns to distinguish polymorphs while maintaining discrete and computationally enumerable inputs. The central methodological benefit of an enumerable input is that it allows users to predict the energy of all possible combinations of spacegroup and Wyckoff positions for a given composition and maximum unit cell size. The lowest-ranked protostructures can then be fed into downstream analysis or modeling. We trained an ensemble of 10 Wrenformer models on the MP-crystals-2022.10.28 dataset for this benchmark.

+ #strong[BOWSR] @zuo_accelerating_2021 (E) - BOWSR combines a symmetry-constrained Bayesian optimizer (BO) with a surrogate energy model to perform an iterative exploration-exploitation-based search of the potential energy landscape. Here we use the pre-trained `Eform_MP_2019` MEGNet model @chen_graph_2019 for the energy model as proposed in the original work. The high sample count needed to explore the PES with BO makes this by far the most expensive model tested.

+ #strong[Voronoi RF] @ward_including_2017 (E) - A random forest trained to map a combination of composition-based Magpie features @ward_generalpurpose_2016 and structure-based relaxation-robust Voronoi tessellation features (effective coordination numbers, structural heterogeneity, local environment properties, …) to DFT formation energies. This fingerprint-based model predates most deep learning for materials but significantly improved over earlier fingerprint-based methods such as the Coulomb matrix @rupp_fast_2012 and partial radial distribution function features @schutt_how_2014. It serves as a baseline model to see how much value the learned featurization of deep learning models can extract from the increasingly large corpus of available training data. We trained Voronoi RF on the MP-crystals-2022.10.28 dataset for this benchmark.

= Data Availability
<sec:data-availability>
The Matbench Discovery training set is the latest Materials Project (MP) @jain_commentary_2013 database release (#link("https://docs.materialsproject.org/changes/database-versions")[v2022.10.28] at time of writing). The test set is the WBM dataset @wang_predicting_2021 which has been archived on Figshare @riebesell_matbenchdata_2023 at #link("https://figshare.com/articles/dataset/22715158"). A snapshot of every ionic step including energies, forces, stresses and magnetic moments in the MP database is also available on Figshare, as are all other data files such as phase diagrams and structures in both ASE and pymatgen format at the same URL #link("https://figshare.com/articles/dataset/22715158").

= Code Availability
<sec:code-availability>
The Matbench Discovery framework, including benchmark implementation, evaluation code, and model submission tools, is available as an open-source Python package and GitHub repository at #link("https://github.com/janosh/matbench-discovery") with a permanent version archived at Zenodo @riebesell_matbenchzenodo_2024. We welcome further model submissions via pull requests.

= Acknowledgements
<sec:acknowledgements>
J.R. acknowledges support from the German Academic Scholarship Foundation (#link("https://wikipedia.org/wiki/Studienstiftung")[Studienstiftung]). A.A.L. acknowledges support from the Royal Society. A.J. and K.A.P. acknowledge the US Department of Energy, Office of Science, Office of Basic Energy Sciences, Materials Sciences and Engineering Division under contract no. DE-AC02-05-CH11231 (Materials Project program KC23MP). This work used computational resources provided by the National Energy Research Scientific Computing Center (NERSC), a U.S. Department of Energy Office of Science User Facility operated under Contract No. DE-AC02-05CH11231.

Our profound gratitude extends to Hai-Chen Wang, Silvana Botti and Miguel A. L. Marques for their valuable contribution in crafting and freely sharing the WBM dataset.

We thank Rickard Armiento, Felix A. Faber and Abhijith S. Parackal for helping develop the evaluation procedures for Wren upon which this work builds. We also thank Rokas Elijosius for assisting in the initial implementation of Wrenformer and Mark Neumann, Luis Barroso-Luque and Yutack Park for submitting compliant models to the leaderboard.

We would like to thank Jason Blake Gibson, Shyue Ping Ong, Chi Chen, Tian Xie, Peichen Zhong and Ekin Dogus Cubuk for helpful discussions.

= Author Contributions
<sec:author-contributions>
Janosh Riebesell: Methodology, Software, Data Curation, Investigation (Training: CGCNN, CGCCN+P, Wrenformer, Voronoi RF), Validation, Formal Analysis, Writing -- original draft. Rhys Goodall: Conceptualization, Software, Validation, Formal Analysis, Writing -- original draft, review & editing. Philipp Benner: Software, Investigation (Training: ALIGNN, MACE), Writing -- original draft. Yuan Chiang: Investigation (Training: MACE), Formal Analysis, Writing -- review & editing. Bowen Deng: Data Curation (MPtrj), Investigation (Training: CHGNet), Writing -- review & editing. Gerbrand Ceder: Supervision, Funding Acquisition. Mark Asta: Supervision, Funding Acquisition. Alpha Lee: Supervision. Anubhav Jain: Supervision. Kristin Persson: Supervision, Writing -- review & editing, Funding Acquisition.

= Competing Interests
<sec:competing-interests>
The authors declare no competing interests.

#pagebreak()
#heading(numbering: none)[Appendix]

#counter(heading).update(0) // reset heading counter
#set heading(numbering: "S1.1") // add S prefix to section counters


= Primer on Crystal Structures
<sec:crystal-primer>
As suggested during review, this section provides a brief introduction to crystal structures for readers unfamiliar with crystallography. We outline why they are foundational to materials science, and how they are represented in the Matbench Discovery dataset.

== Unit Cell
<unit-cell>
The unit cell is the basic building block making up a crystal. It is defined by three lattice vectors ($arrow(a)$, $arrow(b)$, $arrow(c)$) and three angles ($alpha$, $beta$, $gamma$) between these vectors determining the shape and volume of the cell. The entire crystal can be constructed by translating the unit cell along all three lattice vectors.

== Crystal Systems and Symmetry
<crystal-systems-and-symmetry>
Crystals are classified into seven crystal systems (cubic, tetragonal, orthorhombic, hexagonal, trigonal, monoclinic, and triclinic) based on their symmetry. These systems are further divided into 230 space groups that describe all possible sets of symmetry operations (rotations, reflections, and translations) an arrangement of points can have in 3d periodic structures.

== Atomic Positions
<atomic-positions>
Within a unit cell, atomic positions are specified by 3 fractional coordinates (values between 0 and 1) along each lattice vector. These positions, combined with atomic types, fully define the crystal structure. Wyckoff positions denote sets of equivalent atomic sites that transform into each other through the crystal's symmetry operations.

== Structure Representation in Databases
<structure-representation-in-databases>
In materials databases like the Materials Project, crystal structures are typically stored with the following information:

- #strong[Crystal structure]: Cell parameters ($a$, $b$, $c$, $alpha$, $beta$, $gamma$)

- #strong[Atomic positions]: List of element type and fractional coordinates of each atom in the unit cell

- #strong[Site properties]: Optional information about atomic sites (oxidation states, magnetic moments, etc.)

== Energy Quantities
<energy-quantities>
Several energy quantities are used to characterize crystal stability:

- #strong[Formation energy] ($E_(upright("form"))$): Energy released or required when a compound forms from its constituent elements in their lowest-energy pure-element states

- #strong[Energy above hull] ($E_(upright("above hull"))$): Energy relative to the thermodynamic convex hull in a chemical system, indicating stability against decomposition to other competing phases. $E_(upright("above hull")) = 0$, the material is said to be on the hull and thermodynamically stable. If $E_(upright("above hull")) > 0$, the material may still exist with varying degrees of metastability. If $E_(upright("above hull")) > 100 upright(" meV/atom")$, material's are very rarely observed in nature.

- #strong[Decomposition energy]: The energy released during decomposition into more stable phases. Can be negative if a material is lower in any energy than any decomposition products on the known convex hull of its chemical system.

== Matbench Discovery Data Format
<matbench-discovery-data-format>
Matbench Discovery records structures along with DFT-calculated energies and model predictions for 257k materials in the WBM test set. Each entry includes:

- Crystal structure information (lattice, sites, unique material ID)

- Ground truth DFT energies (formation energy, energy above hull)

- Model predictions for comparison

This standardized format enables consistent evaluation of machine learning models for crystal stability prediction and materials discovery.

= Metrics on full test set and for 10k materials predicted most stable
<sec:metrics-for-10k-materials-predicted-most-stable>
As discussed in the first author's PhD thesis @riebesell_machine_2024, unlike @tab:metrics-table-uniq-protos in the main text which evaluates model performance on the subset of unique WBM protostructures, @tab:metrics-table includes the full WBM test set of 257k materials. The 44k additional materials in @tab:metrics-table excluded in @tab:metrics-table-uniq-protos comprise 11175 discarded due to having a matching protostructure in MP plus another 32784 materials which are protostructure duplicates of another WBM material with lower energy. The most noteworthy difference between the two tables is a drop in DAF for all models in @tab:metrics-table compared to @tab:metrics-table-uniq-protos. MACE for example achieves a DAF of 3.5 on the full test set (@tab:metrics-table) compared to 3.85 on the subset of 215.5k materials with unique protostructures (@tab:metrics-table-uniq-protos). This $~$10% increase is largely due to a $~$10% decrease in the fraction of materials below the MP convex hull: 15.3% (32,942 out of 215,488) in @tab:metrics-table-uniq-protos vs 16.7% (42,825 out of 256,963) in the full dataset (@tab:metrics-table). Since DAF is the ratio of the model's precision for stability prediction to the prevalence of stable structures, a lower prevalence results in a higher DAF.

While protostructures are non-trivial to match, thus potentially introducing bias into the test set by trying to deduplicate them, we still opted to feature @tab:metrics-table-uniq-protos in the main text since the removal of overlapping protostructures with MP makes it more closely reflect a model's true ability to extrapolate to out-of-domain materials. Protostructure assignment is based on a combination of Aflow-style prototype labels and the chemical system @hicks_aflow_2021 we use the `get_protostructure_label` implemented in #link("https://github.com/janosh/matbench-discovery/blob/f26f1345bea/matbench_discovery/structure/prototype.py#L105")[`matbench-discovery`] @goodall_aviary_2022. This string encodes the crystal's prototype and chemical system and is invariant to symmetry preserving changes to the structure. Consequently structures with different lattice parameters are flagged as protostructure duplicates as they would both be expected to relax to the same ground state in experimental conditions.

#figure(
  pdf-img("figs/metrics-table-only-compliant.pdf", height: 280pt),
  caption: [
    Metrics on the full test set and for 10k materials predicted most stable.
  ],
)<tab:metrics-table>

A real-world discovery campaign is unlikely to validate all stable predictions from a given model as we did in the main text. Presumably, it will rank model predictions from most to least stable and follow that list as far as time and compute budget permits. Assuming that increases in compute resources will allow average future discovery campaigns to grow in scope, we believe 10k model validations to be a reasonable scope for average campaigns. This is what @tab:metrics-table-first-10k simulates by calculating classification and regression metrics for the 10k test set materials predicted to be most stable by each model. We again show dummy performance in the bottom row. Note that each model is now evaluated on a different slice of the data. However, the bottom row still shows dummy performance across the whole dataset.

#figure(
  pdf-img("figs/metrics-table-first-10k-only-compliant.pdf", height: 300pt),
  caption: [
    Metrics on the full test set and for 10k materials predicted most stable.
  ],
)<tab:metrics-table-first-10k>

= ROC Curves
<sec:roc-curves>
A material is classified as stable if the predicted $E_(upright("above hull"))$ is below the stability threshold. Since all models predict $E_(upright("form"))$, they are insensitive to changes in the threshold $t$. The receiver operating characteristic (ROC) curve for each model is plotted in @fig:roc-models. The diagonal 'No skill' line shows the performance of a dummy model that randomly ranks material stability.

#figure(
  pdf-img("figs/roc-models-only-compliant.pdf"),
  caption: [
    Receiver operating characteristic (ROC) curve for each model. The false positive rate (FPR) on the $x$ axis is the fraction of unstable structures classified as stable. The true positive rate (TPR) on the $y$ axis is the fraction of stable structures classified as stable.
  ],
)<fig:roc-models>

= Parity Plots
<sec:parity-plots>
#figure(
  pdf-img("figs/each-parity-models-4x3-only-compliant.pdf"),
  caption: [
    Parity plots of model-predicted energy distance to the convex hull (based on their formation energy predictions) vs DFT ground truth, color-coded by log density of points. Models are sorted left to right and top to bottom by MAE. For parity plots of formation energy predictions, see @fig:e-form-parity-models.
  ],
)<fig:each-parity-models>

#figure(
  pdf-img("figs/e-form-parity-models-4x3-only-compliant.pdf"),
  caption: [
    Parity plots of model-predicted formation energies vs DFT formation energies, color-coded by log density of points. The models are sorted from left to right and top to bottom by MAE. While similar to the parity plots in @fig:each-parity-models which shows the predicted distance to the convex hull vs DFT ground truth, this figure better visualizes the point density due to formation energy's wider spread. We observe broadly the same failure modes with occasional high DFT energy outliers predicted as near 0 formation energy by the models.
  ],
)<fig:e-form-parity-models>

@fig:each-parity-models shows that all models do well for materials far below the convex hull (left side of the plot). Performance for materials far above the convex hull is more varied with occasional underpredictions of the energy of materials far above the convex hull (right side). All models suffer most in the mode of the distribution at $x = 0$.

Two models stand out as anomalous to the general trends.

Wrenformer is the only model with a large number of severe energy overpredictions at $x = 0$ along the positive $y$ axis. We investigated these failure cases in more detail and found these overpredictions to be dominated by spacegroup 71 with poor representation in the training data. Digging into the MP-crystals-2022.10.28 dataset we see that on post analysis there are a large number of $"A"_2$BC compounds in spacegroup 71 in the v2022.10.28 that are trapped in local minima and thus much higher in energy. Given that Wrenformer nominally predicts the energy the global minima conditioned on the protostructure these are unsuitable for training Wrenformer and their inclusion leads to the errors seen on similar, correctly relaxed structures, in the WBM test set.

The other anomalous model is MACE with several severe underpredictions at $x = 0$ along the negative $y$ axis. We investigated these points for common traits in composition or crystal symmetry but noticed no pattern.

Beyond these MACE outliers visible in the plot, MACE exhibits another rare but reproducible type of failure case, in which the final predicted energy after relaxation is off by several orders of magnitude. The largest 'derailed' prediction was $- 10^22$ eV/atom for `wbm-3-31970` (formula $"H"_2$Ir). In each case, the MACE relaxation exhausted the maximum number of ionic steps set to 500 and caused a volume implosion from initial cell volumes of hundreds to relaxed cell volumes of tens of cubic Angstrom. Using the checkpoint trained on the M3GNet dataset which we received from the MACE authors, this failure mode occurred for several hundred of the 250k test set crystals. Using the checkpoint we trained ourselves on the MPtrj dataset, it affects only 44 test crystals, suggesting that these holes in the MACE PES can perhaps be fully plugged by further increasing the training set or even changing the loss function. Further analysis is ongoing. Since these derailed values are easily identified in practice when actually performing a prospective discovery campaign, we excluded them from the MACE parity plat and all other downstream analyses.

= Hull Distance Box plot
<sec:hull-distance-box-plot>
#figure(
  pdf-img("figs/box-hull-dist-errors-only-compliant.pdf"),
  caption: [
    Box plot of interquartile ranges (IQR) of hull distance errors for each model. The whiskers extend to the 5th and 95th percentiles. The horizontal line inside the box shows the median. BOWSR has the highest median error, while Voronoi RF has the highest IQR.
  ],
)<fig:box-hull-dist-errors>

@fig:box-hull-dist-errors shows a box-plot of the errors in predicting the distance to the convex hull. The box demarks the interquartile range and the whiskers show the 5th and 95th percentiles. BOWSR has the largest median error, while Voronoi RF has the largest IQR. Looking at the models we see that there is a systematic trend to underpredict the distance to the convex hull. This trend can be rationalized give the bias in the Materials Project towards stable and therefore low formation energy structures. In contrast, due to it's construction the WBM test set is distributionally distinct with a higher average formation energy. Therefore we propose that models trained on the Materials Project are likely result in such systematic underpredictions and this affect should be able to be addressed by training on additional data from a similar distribution to WBM i.e. Alexandria @schmidt_improving_2024.

= Classification Histograms using Model-Predicted Energies
<sec:classification-histograms-using-model-predicted-energies>
#figure(
  pdf-img("figs/hist-clf-pred-hull-dist-models-4x3-only-compliant.pdf"),
  caption: [
    Distribution of model-predicted hull distance colored by stability classification. Models are sorted from top to bottom by F1 score. The thickness of the red and yellow bands shows how often models misclassify as a function of how far away from the convex hull they place a material.
  ],
)<fig:hist-clf-pred-hull-dist-models>

@fig:hist-clf-pred-hull-dist-models shows histograms of the model predicted distances to the convex hull colored by the stability classification. These plots allow for practitioners to assess how the accuracy varies for materials depending on how far above or below the hull they are predicted to lie. The most pronounced affect we see from the histograms is that the CGCNN+P histogram is more strongly peaked than CGCNN's giving much better agreement with the actual DFT ground truth distribution of hull distances for the test set. This explains why CGCNN+P performs better as a regressor, but also reveals how it can perform simultaneously worse as a classifier. By moving predictions closer to the stability threshold at 0~eV/atom above the hull, even small errors are significant enough to tip a classification over the threshold.

= Measuring extrapolation performance from WBM batch robustness
<sec:extrapolation-performance>
As a reminder, the WBM test set was generated in 5 successive batches, each step applying another element replacement to an MP source structure or a new stable crystal generated in one of the previous replacement rounds. The likelihood by which one element replaces another is governed by ISCD-mined chemical similarity scores for each pair of elements. Naively, one would expect model performance to degrade with increasing batch count, as repeated substitutions should on average 'diffuse' deeper into uncharted regions of material space, requiring the model to extrapolate more. We observe this effect for some models much more than others.

@fig:rolling-mae-vs-hull-dist-wbm-batches-models shows the rolling MAE as a function of distance to the convex hull for each of the 5 WBM rounds of elemental substitution. These plots show a stronger performance decrease for Wrenformer and Voronoi RF than for FFFs and even force-less GNNs with larger errors like ALIGNN, MEGNet and CGCNN.

#figure(
  pdf-img("figs/tile-rolling-mae-batches-4x3-only-compliant.pdf"),
  caption: [
    Rolling MAE as a function of distance to the convex hull for different models. Most models considered show a predictable decrease in performance from batch 1 to batch 5. The effect is larger for some models than others but batches 4 and 5 consistently incur the highest convex hull distance errors for all models. The FFFs show stronger extrapolative performance, as they show minimal deterioration in performance on later batches that move further away from the original MP training distribution. Wrenformer, by contrast, exhibits a pronounced increase in MAE with batch count. We view these plots as a strong indicator that Matbench Discovery is indeed testing out-of-distribution extrapolation performance as it is Occam's razor explanation for the observed model performance drop with increasing batch count.
  ],
)<fig:rolling-mae-vs-hull-dist-wbm-batches-models>

@fig:rolling-mae-vs-hull-dist-wbm-batches-models shows the rolling MAE for different models split out by calculation batch in the WBM data set. MEGNet and CGCNN both incur a higher rolling MAE than Wrenformer across all 5 batches and across most or all of the hull distance range visible in these plots. However, similar to the FFFs, MEGNet and CGCNN exhibit very little degradation for higher batch counts. The fact that higher errors for later batches are specific to Wrenformer and Voronoi RF suggests that training only on composition and coarse-grained structural features (spacegroup and Wyckoff positions in the case of Wrenformer; coordination numbers, local environment properties, etc. in the case of Voronoi RF) alone is insufficient to learn an extrapolatable map of the PES.

Given its strong performance on batch 1, it is possible that given sufficiently diverse training data, Wrenformer could become similarly accurate to the FFFs across the whole PES landscape at substantially less training and inference cost. However, the loss of predictive forces and stresses may make Wrenformer unattractive for certain applications even then.

= Largest Errors vs.~DFT Hull Distance
<sec:largest-errors-vs.-dft-hull-distance>
Given the large variety of models tested, we asked whether any additional insight into the errors can be gained from looking at how the predictions vary between different models. In @fig:scatter-largest-errors-models-mean-vs-true-hull-dist-all we see two distinct groupings emerge when looking at the 200 structures with the largest errors. This clustering is particularly apparent when points are colored by model disagreement.

#figure(
  pdf-img(
    "figs/scatter-largest-errors-models-mean-vs-true-hull-dist-all-only-compliant.pdf",
    width: 350pt,
  ),
  caption: [
    DFT vs predicted hull distance (average over all models) for the 200 largest error structures. Points are colored by model disagreement as measured by the standard deviation in hull distance predictions from different models. Point scale with the number of atoms in the structures. This plot shows that high-error predictions are biased towards predicting too small hull distances. This is unsurprising considering MP training data mainly consists of low-energy structures. There is a strong color separation between the mostly dark blue low-energy bias predictions and the high-error green or red predictions. Blue means models are in good agreement, i.e. all models are wrong together. Red means large-error predictions with little model agreement, i.e. all models are wrong in different ways. Some of the blue points with large errors yet good agreement among models may be accurate ML predictions for a DFT relaxation gone wrong. The dark blue points also tend to be larger corresponding to larger structures where DFT failures are less surprising. This suggests ML model committees might be used to cheaply screen large databases for DFT errors in a high-throughput manner.
  ],
)<fig:scatter-largest-errors-models-mean-vs-true-hull-dist-all>

= Exploratory Data Analysis
<sec:eda>
To give high-level insights into the MP training and WBM test set used in this work, we include element distributions for structures in both datasets (@fig:element-counts-by-occurrence). To show how frame selection from MP structure relaxation affected relative elemental abundance between MP relaxed structures and the snapshots in MPtrj, @fig:element-counts-ratio-by-occurrence shows element occurrence ratios between MPtrj and MP. Similarly, @fig:mp-vs-mp-trj-vs-wbm-arity-hist shows the elements-per-structure distribution of MP, MPtrj and WBM, normalized by dataset size. The mode of all three datasets is 3, but WBM's share of ternary phases is noticeably more peaked than MP's, which includes small numbers of unary and senary phases. @fig:wbm-energy-hists plots the distributions of formation energies, decomposition energies and convex hull distances (with respect to the convex hull spanned by MP materials only) for the MP and WBM data sets. This plot not only gives insight into the nature of the dataset but also emphasizes the increased difficulty of stability vs formation energy prediction arising from the much narrower distribution of convex hull distances compared to the more spread-out formation energy distribution. For WBM we see that the formation energy distribution exhibits much wider spread than the convex hull distance distribution, spanning almost 10~eV/atom vs less than 1~eV/atom spread in the hull distances. This highlights why stability prediction is a much more challenging task than predicting energy of formation. It requires correctly ranking the subtle energy differences between chemically similar compounds in the same chemical system rather than comparing a single material with the reference energies of its constituent elements. DFT has been shown to significantly benefit from the systematic cancellation of errors between chemically similar systems when trying to identify the lowest-lying polymorph @hautier_accuracy_2012. This beneficial cancellation has yet to be conclusively demonstrated for ML stability predictions. So far, only the lack thereof has been shown in @bartel_critical_2020 where they encountered a much more random error distribution among similar chemistries than simulations from first principles.

Finally, believing MPtrj to be an influential dataset for the near-term continued development of foundational force fields, we plot histograms showing the distributions of target values for energies, forces, stresses and magnetic moments in @fig:mp-trj-hists. The histogram in @fig:mp-trj-n-sites-hist shows the distribution of the number of sites in MPtrj structures. The inset displays the same histogram log-scaled y-axis as well as a cumulative line to show that 90% of MPtrj structures contain fewer than 70 sites. @fig:mp-trj-ptable-hists shows the forces and magnetic moments for the MPTrj dataset projected onto different elements.

#figure(
  [
    #pdf-img("figs/mp-element-counts-by-occurrence-log.pdf")
    #v(-4pt)
    #pdf-img("figs/wbm-element-counts-by-occurrence-log.pdf")
  ],
  caption: [
    The number of structures containing a given element in the MP training set and WBM test set @wang_predicting_2021. The WBM test set in relative terms contains noticeably fewer oxides than MP (and, by extension, MPtrj) with just 11% rather than 53% of structures containing oxygen. Made with pymatviz @riebesell_pymatviz_2022.
  ],
)<fig:element-counts-by-occurrence>

#figure(
  [
    #pdf-img("figs/wbm-mp-ratio-element-counts-by-occurrence-normalized.pdf")
    #v(-4pt)
    #pdf-img("figs/mp-trj-mp-ratio-element-counts-by-occurrence-normalized.pdf")
  ],
  gap: 1em,
  caption: [
    *Top* shows the ratio of elements in the WBM test set to the MP training set. The figure makes clear that WBM explores a distinct chemical space to MP with noble-transition metals, post-transition metals, lanthanides, actinides and metalloids seen more frequently in WBM than MP. Similarly, *Bottom* shows the ratio of elements in the MPtrj dataset to the MP training set. We note a slight overabundance of structures containing hydrogen and halides, indicating that more frames were selected from structures containing these elements which might correlate with the number of ionic steps to find their ground states.
  ],
)<fig:element-counts-ratio-by-occurrence>

#figure(
  pdf-img("figs/mp-vs-mp-trj-vs-wbm-arity-hist.pdf", width: 320pt),
  caption: [
    Distributions of unique elements per structure in MP, MPtrj and the WBM test set. The bar heights are normalized by the total number of structures in each dataset. WBM is dominated by ternary phases making up 75% of the dataset followed by about 13% for quaternaries and 12% for binaries. MP has a more even distribution, in particular with more than double the relative share of quaternary phases and a significant number of quinternaries which are almost absent from WBM. Not shown in this plot for visual clarity are 3% of MP structures containing more than 5 elements (up to 9). We also include MPtrj in this plot to show a slight drop in the relative abundance of quinternaries and higher phases vs MP ground states. This may be due to a poor choice of convergence criteria in early MP relaxation workflows that scaled with the size of the structure (see `EDIFF_PER_ATOM` parameter in `pymatgen` VASP input sets), resulting in unconverged large structures with short relaxation trajectories entering the database. Short relaxations would result in fewer frames of such structures selected for MPtrj. This assumes structures of higher arity correlate with larger structures.
  ],
)<fig:mp-vs-mp-trj-vs-wbm-arity-hist>

#figure(
  grid(
    columns: (1fr, 1fr),
    column-gutter: 9pt,
    row-gutter: 3pt,
    subfigure(
      pdf-img("figs/hist-MP-e-form.pdf"),
      caption: [
        MP formation energy distribution
      ],
      label: <fig:hist-mp-e-form-per-atom>,
    ),
    subfigure(
      pdf-img("figs/hist-MP-hull-dist.pdf"),
      caption: [
        MP decomposition energy distribution
      ],
      label: <fig:hist-mp-hull-dist>,
    ),

    subfigure(
      pdf-img("figs/hist-WBM-e-form.pdf"),
      caption: [
        WBM formation energy distribution
      ],
      label: <fig:hist-wbm-e-form-per-atom>,
    ),
    subfigure(
      pdf-img("figs/hist-WBM-hull-dist.pdf"),
      caption: [
        WBM convex hull distance distribution
      ],
      label: <fig:hist-wbm-hull-dist>,
    ),
  ),
  gap: 18pt,
  caption: [
    Distribution of formation energies, decomposition energy and convex hull distances for the MP training set and the WBM test set. In @fig:hist-mp-e-form-per-atom the bimodality in the MP formation energy distribution is due to the MP anion correction scheme @wang_framework_2021 which significantly lowers some formation energies, especially for oxides. The decomposition energy shown in @fig:hist-mp-hull-dist is calculated as defined in @bartel_critical_2020. We note that multiple materials with the same composition can have negative decomposition energies hence the number with negative decomposition energies is $~$43k compared to $~$35k that are on the hull. Looking at WBM in @fig:hist-wbm-e-form-per-atom and @fig:hist-wbm-hull-dist, we see that the distribution of convex hull distances is tighter due to the change of reference and therefore a more discriminative benchmarking task in terms of goodness-of-fit measures like the coefficient of determination.
  ],
)<fig:wbm-energy-hists>

#figure(
  grid(
    columns: (1fr, 1fr, 1fr),
    column-gutter: 3pt,
    row-gutter: 2em,
    subfigure(
      pdf-img("figs/mp-trj-e-form-hist.pdf"),
      caption: [
        MPTrj formation energy distribution
      ],
      label: <fig:mp-trj-e-form-hist>,
      dy: 10pt,
    ),
    subfigure(
      pdf-img("figs/mp-trj-forces-hist.pdf"),
      caption: [
        MPTrj force distribution
      ],
      label: <fig:mp-trj-forces-hist>,
      dy: 10pt,
    ),
    subfigure(
      pdf-img("figs/mp-trj-stresses-hist.pdf"),
      caption: [
        MPTrj stress distribution
      ],
      label: <fig:mp-trj-stresses-hist>,
      dy: 10pt,
    ),

    subfigure(
      pdf-img("figs/mp-trj-magmoms-hist.pdf"),
      caption: [
        MPTrj magnetic moment distribution
      ],
      label: <fig:mp-trj-magmoms-hist>,
      dy: 10pt,
    ),
    subfigure(
      pdf-img("figs/mp-trj-n-sites-hist.pdf"),
      caption: [
        MPTrj structure size distribution
      ],
      label: <fig:mp-trj-n-sites-hist>,
      dy: 10pt,
    ),
  ),
  gap: 18pt,
  caption: [
    Distribution of energies, forces, stresses, magnetic moments, and number of atoms in MPtrj. Comparing the distribution of formation energies to that of the MP dataset we see that the relative heights of the bimodal peaks is shifted as longer relaxation trajectories are seen on average for materials whose energies are adjusted by the anion correction scheme. The inset in @fig:mp-trj-n-sites-hist uses a log-scale to show the tail of the distribution. The green cumulative line in the inset shows that 82% have less than 50 sites and 97% of structures in MPtrj have less than 100 atoms.
  ],
)<fig:mp-trj-hists>

#figure(
  [
    #pdf-img("figs/mp-trj-magmoms-ptable-hists.pdf")
    #v(-4pt)
    #pdf-img("figs/mp-trj-forces-ptable-hists.pdf")
  ],
  caption: [
    Distribution of magnetic moments and forces for each element MPtrj. This data is used as training targets for all interatomic potentials in this work (only CHGNet uses the absolute value of magnetic moments as targets). The number in the top right corner of each element tile counts the number of target values for that element in all of MPtrj. $y$-axes are log-scaled to reveal the tail of high magnetic moments in some elements. ) reveals rare erroneous data points in MPtrj. For instance, has a single-point calculation with a highly unphysical magnetic moment of 17$mu_(upright("B"))$. For visualization purposes, the $y$-axes are again log-scaled and distributions are truncated at 10~eV/Å. Oxygen has the largest outliers with mean absolute forces of up to 160~eV/Å.
  ],
)<fig:mp-trj-ptable-hists>

#bibliography("references.bib")
