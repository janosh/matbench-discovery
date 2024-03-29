---
Title: Matbench Discovery - An automated leaderboard for benchmarking ML energy models at crystal stability prediction
Website: https://mrs.org/meetings-events/spring-meetings-exhibits/2023-mrs-spring-meeting
Submission date: 2022-10-27
Attendance: indicated in person
Keywords: performance, computing, machine learning, informatics, data, database
Symposium: MD01 Integrating ML with Materials Science
Conference dates: April 10-14 2023
---

The past 2 years have seen the release of several machine learning models specifically designed to predict the DFT final energy given only the unrelaxed crystal structure. Such models are suited for a materials discovery workflow in which they pre-filter and/or pre-relax structures that are then fed into high-throughput DFT. Examples of such models include

- [BOWSR](https://sciencedirect.com/science/article/pii/S1369702121002984): Bayesian optimizer taking crystal symmetries into account when exploring/exploiting known samples of a crystal's potential energy surface
- [Wren](https://www.science.org/doi/10.1126/sciadv.abn4117) [[arxiv](https://arxiv.org/abs/2106.11132)]: uses composition, spacegroup and Wyckoff symmetries as coarse-grained, relaxation-invariant, enumerable descriptors of a crystal structure to predict the post-relaxation energy
- [M3GNet](https://arxiv.org/abs/2202.02450): trained on VASP relaxation trajectories, M3GNet predicts forces and stresses acting within a structure to perform pseudo-relaxation that yields the energy of a structure more closely resembling the DFT-relaxed crystal

However, even closely following the literature, it remains unclear which methodology and in particular which model performs best.
Matbench Discovery extends the growing Matbench ecosystem by building on the dataset published in Wang 2021 npj Computational Materials. This work generated ~250k structures with chemical similarity-based element substitution. After DFT relaxation, they found 10% (~20k) to lie on the Materials Project convex hull. Our benchmark measures how precisely models identify the 20k stable structures from the 250k search pool, as well as how many false positives each one incurs. The latter lead to wasted compute for relaxing unstable structures.

One unusual but compelling feature of the test set is that it was generated by 5 rounds of elemental substitution, each time feeding newly found stable structures back into the pool. This allows for out-of-distribution testing as repeated substitutions on average increase chemical dissimilarity from the training set. Indeed we see a decline in model performance with substitution count which is more pronounced for some models than others. This becomes an important metric that influences model ranking as any prospective discovery campaign for exceptional materials will require models to operate robustly on out-of-distribution crystals. Models that estimate their predictive uncertainty in good correlation with their error receive bonus points.

We also look at how model errors compare across crystal structures and identify classes of materials on which all models perform poorly, indicating a lack of attention from the community and/or lower-quality training data in certain areas of materials space.

Like the original Matbench, Discovery will grow over time as new models are released. To that end, we publish our train and test data along with code to quickly deploy models on them as a `pip`-installable Python package. This allows model authors to add their models to our automated leaderboard by making pull requests to the [Matbench Discovery repo](https://github.com/janosh/matbench-discovery) with their benchmark results and the code used to generate them. The Discovery website is kept up to date through continuous deployment.

The utility of this work is twofold:

1. Make it easy for future researchers to identify SOTA ML models that predict stability given only the unrelaxed crystal. This enables informed decisions when choosing ML filters for high-throughput discovery projects.
2. On a more philosophical note, by keeping the leaderboard up to date with new model releases, we ultimately hope to answer the question of which methodology works best; DFT emulators like M3GNet, one-shot predictors like Wren or a new method yet to be released?
