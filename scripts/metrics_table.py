# %%
from __future__ import annotations

from typing import Any

import pandas as pd
import requests
import wandb
from sklearn.metrics import f1_score, r2_score
from tqdm import tqdm

from matbench_discovery import ROOT, WANDB_PATH, today
from matbench_discovery.data import load_df_wbm_with_preds

__author__ = "Janosh Riebesell"
__date__ = "2022-11-28"


# %%
models: dict[str, dict[str, Any]] = {
    "Wren": dict(n_runs=0),
    "CGCNN": dict(
        n_runs=10,
        filters=dict(
            created_at={"$gt": "2022-11-21", "$lt": "2022-11-23"},
            display_name={"$regex": "cgcnn-robust-formation_energy_per_atom"},
        ),
    ),
    "Voronoi RF": dict(
        n_runs=70,
        filters=dict(
            created_at={"$gt": "2022-11-17", "$lt": "2022-11-28"},
            display_name={"$regex": "voronoi-features"},
        ),
    ),
    "Wrenformer": dict(
        n_runs=10,
        filters=dict(
            created_at={"$gt": "2022-11-14", "$lt": "2022-11-16"},
            display_name={"$regex": "wrenformer-robust-mp-formation_energy"},
        ),
    ),
    "MEGNet": dict(
        n_runs=1,
        filters=dict(
            created_at={"$gt": "2022-11-17", "$lt": "2022-11-19"},
            display_name={"$regex": "megnet-wbm-IS2RE"},
        ),
    ),
    "M3GNet": dict(
        n_runs=99,
        filters=dict(
            created_at={"$gt": "2022-10-31", "$lt": "2022-11-01"},
            display_name={"$regex": "m3gnet-wbm-IS2RE"},
        ),
    ),
    "BOWSR MEGNet": dict(
        n_runs=1000,
        filters=dict(
            created_at={"$gt": "2022-11-22", "$lt": "2022-11-25"},
            display_name={"$regex": "bowsr-megnet"},
        ),
    ),
}

run_times: dict[str, dict[str, str | int | float]] = {}


# %% calculate total model run times from wandb logs
# NOTE these model run times are pretty meaningless since some models were run on GPU
# (Wrenformer and CGCNN), others on CPU. Also BOWSR MEGNet, M3GNet and MEGNet weren't
# trained from scratch. Their run times only indicate the time needed to predict the
# test set.

for model in (pbar := tqdm(models)):
    model_dict = models[model]
    n_runs, filters = (model_dict.get(x) for x in ("n_runs", "filters"))
    if n_runs == 0 or model in run_times:
        continue
    pbar.set_description(model)

    runs = wandb.Api().runs(WANDB_PATH, filters=filters)

    assert len(runs) == n_runs, f"found {len(runs)=} for {model}, expected {n_runs}"

    run_time = sum(run.summary.get("_wandb", {}).get("runtime", 0) for run in runs)
    # NOTE we assume all jobs have the same metadata here
    metadata = requests.get(runs[0].file("wandb-metadata.json").url).json()

    n_gpu, n_cpu = metadata.get("gpu_count", 0), metadata.get("cpu_count", 0)
    run_times[model] = {"Run time": run_time, "Hardware": f"GPU: {n_gpu}, CPU: {n_cpu}"}


# on 2022-11-28:
# run_times = {'Voronoi RF': 739608,
#  'Wrenformer': 208399,
#  'MEGNet': 12396,
#  'M3GNet': 301138,
#  'BOWSR MEGNet': 9105237}


# %%
df_wbm: pd.DataFrame = load_df_wbm_with_preds(models=list(models))
target_col = "e_form_per_atom_mp2020_corrected"
df_wbm = df_wbm.round(3).query(f"{target_col} < 5")

e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"
e_above_hull = df_wbm[e_above_hull_col]


# %%
df_metrics = pd.DataFrame(run_times).T

for model in models:
    dct = {}
    e_above_hull_pred = df_wbm[model] - df_wbm[target_col]

    dct["F1"] = f1_score(e_above_hull < 0, e_above_hull_pred < 0)
    dct["Precision"] = f1_score(e_above_hull < 0, e_above_hull_pred < 0, pos_label=True)
    dct["Recall"] = f1_score(e_above_hull < 0, e_above_hull_pred < 0, pos_label=False)

    dct["MAE"] = (e_above_hull_pred - e_above_hull).abs().mean()

    dct["RMSE"] = ((e_above_hull_pred - e_above_hull) ** 2).mean() ** 0.5
    dct["R2"] = r2_score(
        e_above_hull.loc[e_above_hull_pred.dropna().index], e_above_hull_pred.dropna()
    )

    df_metrics.loc[model, list(dct)] = dct.values()


df_styled = df_metrics.style.format(precision=3).background_gradient(
    cmap="viridis",
    # gmap=np.log10(df_table) # for log scaled color map
)


# %%
styles = {
    "": "font-family: sans-serif; border-collapse: collapse;",
    "td, th": "border: 1px solid #ddd; text-align: left; padding: 8px;",
}
df_styled.set_table_styles([dict(selector=sel, props=styles[sel]) for sel in styles])

html_path = f"{ROOT}/figures/{today}-metrics-table.html"
# df_styled.to_html(html_path)
