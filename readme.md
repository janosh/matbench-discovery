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

> TL;DR: We benchmark ML models on crystal stability prediction from unrelaxed structures finding interatomic potentials in particular to be a valuable addition to high-throughput discovery pipelines.

Matbench Discovery is an [interactive leaderboard](https://janosh.github.io/matbench-discovery) and associated [PyPI package](https://pypi.org/project/matbench-discovery) which together make it easy to benchmark ML energy models on a task designed to closely simulate a high-throughput discovery campaign for new stable inorganic crystals.

In version 1 of this benchmark, we explore 8 models covering multiple methodologies ranging from random forests to graph neural networks, from one-shot predictors to iterative Bayesian optimizers and interatomic potential-based relaxers. We find [CHGNet](https://github.com/CederGroupHub/chgnet) ([paper](https://doi.org/10.48550/arXiv.2302.14231)) to achieve the highest F1 score of 0.59, $R^2$ of 0.61 and a discovery acceleration factor (DAF) of 3.06 (meaning a 3x higher rate of stable structures compared to dummy selection in our already enriched search space). See the [**full results**](https://janosh.github.io/matbench-discovery/paper#results) in our interactive dashboard which provides valuable insights for maintainers of large-scale materials databases. We show these models have become powerful enough to warrant deploying them as triaging steps to more effectively allocate compute in high-throughput DFT relaxations.

<slot name="metrics-table" />

We welcome contributions that add new models to the leaderboard through [GitHub PRs](https://github.com/janosh/matbench-discovery/pulls). See the [usage and contributing guide](https://janosh.github.io/matbench-discovery/contribute) for details.

For a version 2 release of this benchmark, we plan to merge the current training and test sets into the new training set and acquire a much larger test set (potentially at meta-GGA level of theory) compared to the v1 test set of 257k structures. Anyone interested in joining this effort please [open a GitHub discussion](https://github.com/janosh/matbench-discovery/discussions) or [reach out privately](mailto:janosh@lbl.gov?subject=Matbench%20Discovery).

For detailed results and analysis, check out the [paper](https://janosh.github.io/matbench-discovery/paper) and [supplementary material](https://janosh.github.io/matbench-discovery/si).
