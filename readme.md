<slot name="title">
  <h1 align="center">
    <img src="https://github.com/janosh/matbench-discovery/raw/main/site/static/favicon.svg" alt="Logo" width="60px"><br>
    Matbench Discovery
  </h1>
</slot>

<h4 align="center" class="toc-exclude" style="display: none;">

[![arXiv](https://img.shields.io/badge/arXiv-2308.14920-blue?logo=arxiv&logoColor=white)](https://arxiv.org/abs/2308.14920)
[![Tests](https://github.com/janosh/matbench-discovery/actions/workflows/test.yml/badge.svg)](https://github.com/janosh/matbench-discovery/actions/workflows/test.yml)
[![GitHub Pages](https://github.com/janosh/matbench-discovery/actions/workflows/gh-pages.yml/badge.svg)](https://github.com/janosh/matbench-discovery/actions/workflows/gh-pages.yml)
[![Requires Python 3.11+](https://img.shields.io/badge/Python-3.11+-blue.svg?logo=python&logoColor=white)](https://python.org/downloads)
[![PyPI](https://img.shields.io/pypi/v/matbench-discovery?logo=pypi&logoColor=white)](https://pypi.org/project/matbench-discovery?logo=pypi&logoColor=white)

</h4>

Matbench Discovery is an [interactive leaderboard](https://janosh.github.io/matbench-discovery/models) which ranks ML models on multiple tasks designed to simulate high-throughput discovery of new stable inorganic crystals, finding their ground state atomic positions and predicting their thermal conductivity.

We rank <slot name="model_count">20+</slot> models covering multiple methodologies including graph neural network (GNN) interatomic potentials, GNN one-shot predictors, iterative Bayesian optimizers and random forests with shallow-learning structure fingerprints.

<slot name="best_report" />

Our results show that ML models have become robust enough to deploy them as triaging steps to more effectively allocate compute in high-throughput DFT relaxations. This work provides valuable insights for anyone looking to build large-scale materials databases.

To cite Matbench Discovery, use:

> Janosh Riebesell, Rhys E. A. Goodall, Philipp Benner, Yuan Chiang, Bowen Deng, Mark Asta, Gerbrand Ceder, Alpha A. Lee, Anubhav Jain, and Kristin A. Persson. "Matbench Discovery -- A Framework to Evaluate Machine Learning Crystal Stability Predictions." arXiv, August 28, 2023. https://doi.org/10.1038/s42256-025-01055-1.

We welcome new models additions to the leaderboard through GitHub PRs. See the [contributing guide](https://janosh.github.io/matbench-discovery/contribute) for details and ask support questions via [GitHub discussion](https://github.com/janosh/matbench-discovery/discussions).

For detailed results and analysis, check out https://nature.com/articles/s42256-025-01055-1.

> Disclaimer: We evaluate how accurately ML models predict several material properties like thermodynamic stability, thermal conductivity, and atomic positions, in all cases using PBE DFT as reference data. Although these properties are important for high-throughput materials discovery, the ranking cannot give a complete picture of a model's overall ability to drive materials research. A high ranking does not constitute endorsement by the Materials Project.
