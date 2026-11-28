The [`matbench-discovery`](https://pypi.org/project/matbench-discovery) Python package is this benchmark's reference implementation: it fetches the training and test sets, provides the shared task runners every model submission plugs into, and computes each leaderboard metric.

Install into an activated Python 3.14+ virtual environment:

```sh
uv pip install matbench-discovery
```

## Quickstart

Load the WBM test set with a model's predictions and reproduce its discovery metrics:

```python
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import Model, TestSubset
from matbench_discovery.metrics import discovery

model = Model.mace_mp_0

# downloads (and caches) the WBM test set plus this model's predictions
df_wbm = load_df_wbm_with_preds(models=[model])

# mask outliers and round like the leaderboard, then score against the fixed DFT hull
df_ref, preds = discovery.prepare_model_predictions(df_wbm, df_wbm[model.key])
metrics = discovery.calc_discovery_metrics(
    df_ref, preds, uniq_proto_prevalence=discovery.wbm_uniq_proto_prevalence()
)
print(metrics[TestSubset.full_test_set])
# {'F1': 0.668, 'DAF': 3.411, 'Precision': 0.584, 'Recall': 0.78, ...}
```

Downloads are cached under `~/.cache/matbench-discovery` (override with `MBD_CACHE_DIR`).

## Where to look

- **Registries and data**: `enums` holds the `Model` and `DataFiles` registries plus the column-name enums (`Key`, `Task`, `TestSubset`); `data` loads them into dataframes.
- **Metrics**: `metrics.discovery`, `metrics.geo_opt`, `metrics.phonons`, `metrics.md` and `metrics.diatomics` compute each task's leaderboard numbers.
- **Task runners**: `discovery`, `phonons`, `md` and `diatomics` drive the evaluations, `calculators` builds each model's ASE calculator, and `ase_relax` performs the relaxations.
- **Supporting modules**: `figs` and `data_figs` build this site's plots, `remote` handles downloads and Figshare uploads, `cli` and `runner_cli` parse the eval-script and shared-runner flags, and `hpc` handles Slurm submission, sharding and run-cost provenance.

Contributing a model means registering a calculator, not calling these modules yourself (see [the contribution guide](/contribute)).

Everything below is generated from the package's docstrings, and each `source` badge links to the exact commit this site was built from.
