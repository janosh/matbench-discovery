<h1 align="center">
  <img src="https://github.com/janosh/matbench-discovery/raw/main/site/static/favicon.svg" alt="Logo" width="60px"><br>
  Matbench Discovery
</h1>

<h4 align="center" class="toc-exclude">

[![arXiv](https://img.shields.io/badge/arXiv-2308.14920-blue?logo=arxiv&logoColor=white)](https://arxiv.org/abs/2308.14920)
[![Tests](https://github.com/janosh/matbench-discovery/actions/workflows/test.yml/badge.svg)](https://github.com/janosh/matbench-discovery/actions/workflows/test.yml)
[![GitHub Pages](https://github.com/janosh/matbench-discovery/actions/workflows/gh-pages.yml/badge.svg)](https://github.com/janosh/matbench-discovery/actions/workflows/gh-pages.yml)
[![Requires Python 3.11+](https://img.shields.io/badge/Python-3.11+-blue.svg?logo=python&logoColor=white)](https://python.org/downloads)
[![PyPI](https://img.shields.io/pypi/v/matbench-discovery?logo=pypi&logoColor=white)](https://pypi.org/project/matbench-discovery?logo=pypi&logoColor=white)

</h4>

> TL;DR: We benchmark ML models on crystal stability prediction from unrelaxed structures finding universal interatomic potentials (UIP) like [CHGNet](https://github.com/CederGroupHub/chgnet), [MACE](https://github.com/ACEsuit/mace) and [M3GNet](https://github.com/materialsvirtuallab/m3gnet) to be highly accurate, robust across chemistries and ready for production use in high-throughput materials discovery.

Matbench Discovery is an [interactive leaderboard](https://janosh.github.io/matbench-discovery/models) and associated [PyPI package](https://pypi.org/project/matbench-discovery) which together make it easy to rank ML energy models on a task designed to simulate a high-throughput discovery campaign for new stable inorganic crystals.

We've tested <slot name="model-count" />models covering multiple methodologies ranging from random forests with structure fingerprints to graph neural networks, from one-shot predictors to iterative Bayesian optimizers and interatomic potential relaxers.

<slot name="best-report" />

Our results show that ML models have become robust enough to deploy them as triaging steps to more effectively allocate compute in high-throughput DFT relaxations. This work provides valuable insights for anyone looking to build large-scale materials databases.

<slot name="metrics-table" />

We welcome contributions that add new models to the leaderboard through GitHub PRs. See the [contributing guide](https://janosh.github.io/matbench-discovery/contribute) for details.

If you're interested in joining this work, feel free to [open a GitHub discussion](https://github.com/janosh/matbench-discovery/discussions) or [send an email](mailto:janosh.riebesell@gmail.gov?subject=Collaborate%20on%20Matbench%20Discovery).

For detailed results and analysis, check out the [preprint](https://janosh.github.io/matbench-discovery/preprint).
