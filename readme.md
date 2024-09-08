<h1 align="center" style="line-height: 0; margin: 0 auto 1em;">
  <img src="https://github.com/janosh/matbench-discovery/raw/main/site/static/favicon.svg" alt="Logo" width="60px"><br>
  Matbench Discovery
</h1>

<h4 align="center" class="toc-exclude" style="display: none;">

[![arXiv](https://img.shields.io/badge/arXiv-2308.14920-blue?logo=arxiv&logoColor=white)](https://arxiv.org/abs/2308.14920)
[![Tests](https://github.com/janosh/matbench-discovery/actions/workflows/test.yml/badge.svg)](https://github.com/janosh/matbench-discovery/actions/workflows/test.yml)
[![GitHub Pages](https://github.com/janosh/matbench-discovery/actions/workflows/gh-pages.yml/badge.svg)](https://github.com/janosh/matbench-discovery/actions/workflows/gh-pages.yml)
[![Requires Python 3.11+](https://img.shields.io/badge/Python-3.11+-blue.svg?logo=python&logoColor=white)](https://python.org/downloads)
[![PyPI](https://img.shields.io/pypi/v/matbench-discovery?logo=pypi&logoColor=white)](https://pypi.org/project/matbench-discovery?logo=pypi&logoColor=white)

</h4>

> Disclaimer: We evaluate how accurately ML models predict solid-state thermodynamic stability. Although this is an important aspect of high-throughput materials discovery, the ranking cannot give a complete picture of a model's general applicability to materials. A high ranking does not constitute endorsement by the Materials Project.

Matbench Discovery is an [interactive leaderboard](https://janosh.github.io/matbench-discovery/models) and associated [PyPI package](https://pypi.org/project/matbench-discovery) which together make it easy to rank ML energy models on a task designed to simulate high-throughput discovery of new stable inorganic crystals.

We've tested <slot name="model-count" />models covering multiple methodologies including graph neural network (GNN) interatomic potentials, GNN one-shot predictors, iterative Bayesian optimizers and random forests with shallow-learning structure fingerprints.

<slot name="best-report" />

Our results show that ML models have become robust enough to deploy them as triaging steps to more effectively allocate compute in high-throughput DFT relaxations. This work provides valuable insights for anyone looking to build large-scale materials databases.

<slot name="metrics-table" />

If you'd like to refer to Matbench Discovery in a publication, please cite the [preprint](https://doi.org/10.48550/arXiv.2308.14920):

> Janosh Riebesell, Rhys E. A. Goodall, Philipp Benner, Yuan Chiang, Bowen Deng, Alpha A. Lee, Anubhav Jain, and Kristin A. Persson. "Matbench Discovery -- A Framework to Evaluate Machine Learning Crystal Stability Predictions." arXiv, August 28, 2023. https://doi.org/10.48550/arXiv.2308.14920.

We welcome new models additions to the leaderboard through GitHub PRs. See the [contributing guide](https://janosh.github.io/matbench-discovery/contribute) for details and ask support questions via [GitHub discussion](https://github.com/janosh/matbench-discovery/discussions).

For detailed results and analysis, check out the [preprint](https://janosh.github.io/matbench-discovery/preprint).
