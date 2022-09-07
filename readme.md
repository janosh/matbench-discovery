# Matbench Discovery

[![Link check](https://github.com/janosh/matbench-discovery/actions/workflows/test.yml/badge.svg)](https://github.com/janosh/matbench-discovery/actions/workflows/test.yml)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/janosh/matbench-discovery/main.svg?badge_token=Qza33izjRxSbegTqeSyDvA)](https://results.pre-commit.ci/latest/github/janosh/matbench-discovery/main?badge_token=Qza33izjRxSbegTqeSyDvA)
[![Requires Python 3.9+](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://python.org/downloads)

Several new energy models specifically designed to handle unrelaxed structures were published in 2021/22

- [BOWSR](https://sciencedirect.com/science/article/pii/S1369702121002984)
- [M3GNet](https://arxiv.org/abs/2202.02450)
- [Wren](https://arxiv.org/abs/2106.11132)
- missing one? Please open an issue.

Such models are suited for a materials discovery workflow in which they pre-filter and/or pre-relax structures that are then fed into high-throughput DFT. Even for someone trying to keep up with the literature though, it's unclear which model performs best at that task. Consequently, we think a follow up paper to the 2020 work from Chris Bartel is in order.

[A critical examination of compound stability predictions from machine-learned formation energies](https://nature.com/articles/s41524-020-00362-y)

This project aims to complement Matbench using the WBM dataset published in [Predicting stable crystalline compounds using chemical similarity](https://nature.com/articles/s41524-020-00481-6). They generated ~250k structures with chemical similarity-based elemental substitution and relaxed all of them. ~20k or 10% were found to lie on the Materials Project convex hull. They did 5 iterations of this substitution process. This is a unique and compelling feature of the dataset as it allows out-of distribution testing. We can look at how a model performs when asked to predict on structures increasingly more different from the training set (which is restricted to MP for all models in this benchmark at the moment) since repeated substitutions should - on average - increase chemical dissimilarity.

A good set of baseline models would be CGCNN, Wren and Voronoi tessellation combined with a random forest.
In addition to CGCNN, Wren and Voronoi plus RF, this benchmark includes ALIGNN (current Matbench SOTA), BOWSR and M3GNet to see how many of the 20k stable structures each of these models recovers and how their performance changes as a function of iteration number, i.e. how well they extrapolate. Like Matbench, future model submission to this benchmark can be added via PRs to this repo.

Soon to come: The performance of all models will be displayed in an online dashboard similar to <https://matbench.materialsproject.org>. Our goal is to make it easy for future research to compare new models to the SOTA and help potential users make an informed decisions which model to go with.

On a more philosophical note: Another primary goal of this benchmark is to at least partly answer the question how useful ML energy models really are at helping to accelerate inorganic crystal searching and whether DFT emulators like M3GNet or one-shot predictors like Wren do better.
