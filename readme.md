<h1 align="center" style="display: grid;">
<img src="https://raw.githubusercontent.com/janosh/matbench-discovery/main/site/static/favicon.svg" alt="Logo" width="80px">
Matbench Discovery
</h1>

<h4 align="center" class="toc-exclude">

[![Tests](https://github.com/janosh/matbench-discovery/actions/workflows/test.yml/badge.svg)](https://github.com/janosh/matbench-discovery/actions/workflows/test.yml)
[![GitHub Pages](https://github.com/janosh/matbench-discovery/actions/workflows/gh-pages.yml/badge.svg)](https://github.com/janosh/matbench-discovery/actions/workflows/gh-pages.yml)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/janosh/matbench-discovery/main.svg?badge_token=Qza33izjRxSbegTqeSyDvA)](https://results.pre-commit.ci/latest/github/janosh/matbench-discovery/main?badge_token=Qza33izjRxSbegTqeSyDvA)
[![Requires Python 3.9+](https://img.shields.io/badge/Python-3.9+-blue.svg?logo=python&logoColor=white)](https://python.org/downloads)
[![PyPI](https://img.shields.io/pypi/v/matbench-discovery?logo=pypi&logoColor=white)](https://pypi.org/project/matbench-discovery?logo=pypi&logoColor=white)

</h4>

Matbench Discovery is an [interactive leaderboard](https://janosh.github.io/matbench-discovery) and associated [PyPI package](https://pypi.org/project/matbench-discovery) for benchmarking ML energy models on a task designed to closely emulate a real-world computational materials discovery workflow. In it, these models take on the role of a triaging step prior to DFT to decide how to allocate limited compute budget for structure relaxations.

We welcome contributions that add new models to the leaderboard through [GitHub PRs](https://github.com/janosh/matbench-discovery/pulls). See the [usage and contributing guide](https://janosh.github.io/matbench-discovery/how-to-contribute) for details.

Several new energy models specifically designed to handle unrelaxed structures were published in 2021/22

- [BOWSR](https://sciencedirect.com/science/article/pii/S1369702121002984)
- [M3GNet](https://arxiv.org/abs/2202.02450)
- [Wren](https://arxiv.org/abs/2106.11132)
- missing one? Please open an issue.

Such models are suited for a materials discovery workflow in which they pre-filter and/or pre-relax structures that are then fed into high-throughput DFT. Even for someone trying to keep up with the literature though, it's unclear which model performs best at that task. Consequently, we think a follow-up paper to the 2020 work from Chris Bartel is in order.

[A critical examination of compound stability predictions from machine-learned formation energies](https://nature.com/articles/s41524-020-00362-y)

This project aims to complement Matbench using the **WBM dataset** published in [Predicting stable crystalline compounds using chemical similarity](https://nature.com/articles/s41524-020-00481-6). They generated ~250k structures with chemical similarity-based elemental substitution and relaxed all of them. ~20k or 10% were found to lie on the Materials Project convex hull. They did 5 iterations of this substitution process. This is a unique and compelling feature of the dataset as it allows out-of-distribution testing. We can look at how a model performs when asked to predict on structures increasingly more different from the training set (which is restricted to MP for all models in this benchmark at the moment) since repeated substitutions should - on average - increase chemical dissimilarity.

A good set of baseline models would be CGCNN, Wren and Voronoi tessellation combined with a random forest. In addition to CGCNN, Wren and Voronoi plus RF, this benchmark includes BOWSR and M3GNet to see how many of the 20k stable structures each of these models recover and how their performance changes as a function of iteration number, i.e. how well they extrapolate. Like Matbench, future model submissions to this benchmark can be added via PRs to this repo.

Our goal with this site is to serve as an interactive dashboard for researchers that makes it easy to compare the performance of different energy models on metrics like precision, recall and discovery acceleration to find the model that best suits your needs. You can then make an informed decision about which model to pick by trading off compute savings from an increased hit rate to a more complete discovery in your materials space of interest from higher recall.

On a more philosophical note: Another primary goal of this benchmark is to at least partly answer the question of how useful ML energy models really are at helping to accelerate inorganic crystal searching and whether DFT emulators like M3GNet or one-shot predictors like Wren do better.
