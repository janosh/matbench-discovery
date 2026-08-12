The [`matbench-discovery`](https://pypi.org/project/matbench-discovery) Python package is the reference implementation behind this benchmark. It fetches the training and test sets, provides the shared task runners that every model submission plugs into, and computes each metric shown on the leaderboards.

```sh
pip install matbench-discovery
```

## Quickstart

Load the WBM test set alongside a model's predictions and reproduce its discovery metrics:

```python
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.metrics.discovery import stable_metrics

model = Model.mace_mp_0

# downloads (and caches) the WBM test set plus this model's predictions
df_wbm = load_df_wbm_with_preds(models=[model])

# predicted hull distance = DFT hull distance + the model's formation energy error
each_pred = df_wbm[MbdKey.each_true] + df_wbm[model.label] - df_wbm[MbdKey.e_form_dft]

print(stable_metrics(df_wbm[MbdKey.each_true], each_pred))
# {'F1': 0.665, 'DAF': 3.457, 'Precision': 0.576, 'Recall': 0.787, ...}
```

Downloads are cached in `data/` if you cloned the repo and under `~/.cache/matbench-discovery` otherwise. Set `MBD_CACHE_DIR` to override.

## Where to look

- **Registries and data** — `enums` holds the `Model` and `DataFiles` registries plus the column-name enums (`MbdKey`, `Task`, `TestSubset`), and `data` loads them into dataframes.
- **Metrics** — `metrics.discovery`, `metrics.geo_opt`, `metrics.phonons`, `metrics.md` and `metrics.diatomics` compute the leaderboard numbers for each task.
- **Task runners** — `discovery`, `phonons`, `md` and `diatomics` drive the evaluations, `calculators` builds each model's ASE calculator, and `ase_relax` performs the relaxations.
- **Supporting modules** — `figs` and `data_figs` build this site's plots, `remote` handles downloads and Figshare uploads, and `hpc`, `cli` and `runner_cli` back the Slurm submission scripts.

Contributing a model means registering a calculator rather than calling these modules yourself — see [the contribution guide](/contribute).

Everything below is generated from the package's docstrings, and each `source` badge links to the exact commit this site was built from.
