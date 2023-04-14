<h1 align="center">
  <img src="https://github.com/janosh/matbench-discovery/raw/main/site/static/favicon.svg" alt="Logo" width="60px"><br>
  Matbench Discovery
</h1>

<h4 align="center" class="toc-exclude">

[![Tests](https://github.com/janosh/matbench-discovery/actions/workflows/test.yml/badge.svg)](https://github.com/janosh/matbench-discovery/actions/workflows/test.yml)
[![GitHub Pages](https://github.com/janosh/matbench-discovery/actions/workflows/gh-pages.yml/badge.svg)](https://github.com/janosh/matbench-discovery/actions/workflows/gh-pages.yml)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/janosh/matbench-discovery/main.svg?badge_token=Qza33izjRxSbegTqeSyDvA)](https://results.pre-commit.ci/latest/github/janosh/matbench-discovery/main?badge_token=Qza33izjRxSbegTqeSyDvA)
[![Requires Python 3.9+](https://img.shields.io/badge/Python-3.9+-blue.svg?logo=python&logoColor=white)](https://python.org/downloads)
[![PyPI](https://img.shields.io/pypi/v/matbench-discovery?logo=pypi&logoColor=white)](https://pypi.org/project/matbench-discovery?logo=pypi&logoColor=white)

</h4>

> TL;DR: We benchmark ML models on crystal stability prediction from unrelaxed structures finding interatomic potentials in particular to be a valuable addition to high-throughput discovery pipelines.

Matbench Discovery is an [interactive leaderboard](https://janosh.github.io/matbench-discovery/models) and associated [PyPI package](https://pypi.org/project/matbench-discovery) which together make it easy to benchmark ML energy models on a task designed to closely simulate a high-throughput discovery campaign for new stable inorganic crystals.

So far, we've tested 8 models covering multiple methodologies ranging from random forests with structure fingerprints to graph neural networks, from one-shot predictors to iterative Bayesian optimizers and interatomic potential-based relaxers. We find [CHGNet](https://github.com/CederGroupHub/chgnet) ([paper](https://doi.org/10.48550/arXiv.2302.14231)) to achieve the highest F1 score of 0.59, $R^2$ of 0.61 and a discovery acceleration factor (DAF) of 3.06 (meaning a 3x higher rate of stable structures compared to dummy selection in our already enriched search space). We believe our results show that ML models have become robust enough to deploy them as triaging steps to more effectively allocate compute in high-throughput DFT relaxations. This work provides valuable insights for anyone looking to build large-scale materials databases.

<slot name="metrics-table" />

We welcome contributions that add new models to the leaderboard through GitHub PRs. See the [contributing guide](https://janosh.github.io/matbench-discovery/contribute) for details.

Anyone interested in joining this effort please [open a GitHub discussion](https://github.com/janosh/matbench-discovery/discussions) or [reach out privately](mailto:janosh@lbl.gov?subject=Matbench%20Discovery).

For detailed results and analysis, check out our [preprint](https://janosh.github.io/matbench-discovery/preprint) and [SI](https://janosh.github.io/matbench-discovery/si).
