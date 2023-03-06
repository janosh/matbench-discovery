"""Compile metrics and total run times for all models and export them to JSON, a
pandas-styled HTML table and a plotly figure.
"""


# %%
from __future__ import annotations

from typing import Any

import dataframe_image as dfi
import pandas as pd
import requests
import wandb
import wandb.apis.public
from pymatviz.utils import save_fig
from tqdm import tqdm

from matbench_discovery import FIGS, MODELS, ROOT, WANDB_PATH
from matbench_discovery.plots import px
from matbench_discovery.preds import PRED_FILES, df_metrics, df_wbm

__author__ = "Janosh Riebesell"
__date__ = "2022-11-28"


# %%
model_stats: dict[str, dict[str, str | int | float]] = {}
models: dict[str, dict[str, Any]] = {
    "BOWSR + MEGNet": dict(
        n_runs=500,
        filters=dict(
            display_name={"$regex": "bowsr-megnet"},
            created_at={"$gt": "2023-01-20", "$lt": "2023-01-22"},
        ),
    ),
    "CHGNet": dict(
        n_runs=102,
        filters=dict(
            display_name={"$regex": "chgnet-wbm-IS2RE-"},
            created_at={"$gt": "2023-03-05", "$lt": "2023-03-07"},
        ),
    ),
    "CGCNN": dict(
        n_runs=10,
        filters=dict(
            display_name={"$regex": "cgcnn-robust-formation_energy_per_atom"},
            created_at={"$gt": "2022-11-21", "$lt": "2022-11-23"},
        ),
    ),
    "M3GNet": dict(
        n_runs=99,
        filters=dict(
            display_name={"$regex": "m3gnet-wbm-IS2RE"},
            created_at={"$gt": "2022-10-31", "$lt": "2022-11-01"},
        ),
    ),
    "MEGNet": dict(
        n_runs=1,
        filters=dict(
            display_name={"$regex": "megnet-wbm-IS2RE"},
            created_at={"$gt": "2022-11-17", "$lt": "2022-11-19"},
        ),
    ),
    "Voronoi RF": dict(
        n_runs=68,
        filters=dict(
            display_name={"$regex": "voronoi-features"},
            created_at={"$gt": "2022-11-17", "$lt": "2022-11-28"},
        ),
    ),
    "Wrenformer": dict(
        n_runs=10,
        filters=dict(
            display_name={"$regex": "wrenformer-robust-mp-formation_energy"},
            created_at={"$gt": "2022-11-14", "$lt": "2022-11-16"},
        ),
    ),
}

assert not (
    unknown_models := set(models) - set(PRED_FILES)
), f"{unknown_models=} missing predictions file"


# %% calculate total model run times from wandb logs
# NOTE these model run times are pretty meaningless since some models were run on GPU
# (Wrenformer and CGCNN), others on CPU. Also BOWSR + MEGNet, M3GNet and MEGNet weren't
# trained from scratch. Their run times only indicate the time needed to predict the
# test set.

for model in (pbar := tqdm(models)):
    n_runs, filters = (models[model].get(x) for x in ("n_runs", "filters"))
    if n_runs == 0 or model in model_stats:
        continue
    pbar.set_description(model)
    if "runs" in models[model]:
        runs: wandb.apis.public.Runs = models[model]["runs"]
    else:
        models[model]["runs"] = runs = wandb.Api().runs(WANDB_PATH, filters=filters)

    assert len(runs) == n_runs, f"found {len(runs)=} for {model}, expected {n_runs}"

    each_run_time = [run.summary.get("_wandb", {}).get("runtime", 0) for run in runs]

    run_time_total = sum(each_run_time)
    # NOTE we assume all jobs have the same metadata here
    metadata = requests.get(runs[0].file("wandb-metadata.json").url).json()

    n_gpu, n_cpu = metadata.get("gpu_count", 0), metadata.get("cpu_count", 0)
    model_stats[model] = {
        (time_col := "Run Time (h)"): run_time_total / 3600,
        "GPU": n_gpu,
        "CPU": n_cpu,
        "Slurm Jobs": n_runs,
    }


ax = (pd.Series(each_run_time) / 3600).hist(bins=100)
ax.set(
    title=f"Run time distribution for {model}", xlabel="Run time [h]", ylabel="Count"
)

model_stats["M3GNet + MEGNet"] = model_stats["M3GNet"].copy()
model_stats["M3GNet + MEGNet"][time_col] = (
    model_stats["MEGNet"][time_col] + model_stats["M3GNet"][time_col]  # type: ignore
)
model_stats["CGCNN+P"] = {}


df_stats = pd.concat([df_metrics, pd.DataFrame(model_stats)]).T


# %%
higher_is_better = ["DAF", "R²", "Precision", "Recall", "F1", "Accuracy", "TPR", "TNR"]
lower_is_better = ["MAE", "RMSE", "FNR", "FPR"]
styler = (
    df_metrics.T.rename(columns={"R2": "R²"})
    # append arrow up/down to table headers to indicate higher/lower metric is better
    # .rename(columns=lambda x: x + " ↑" if x in higher_is_better else x + " ↓")
    .style.format(precision=2)
    # reverse color map if lower=better
    .background_gradient(cmap="viridis_r", subset=lower_is_better)
    # .background_gradient(
    #     cmap="viridis_r",
    #     subset=[time_col],
    #     gmap=np.log10(df_stats[time_col].to_numpy()),  # for log scaled color map
    # )
    .background_gradient(cmap="viridis", subset=higher_is_better)
)
styles = {
    "": "font-family: sans-serif; border-collapse: collapse;",
    "td, th": "border: none; padding: 4px 6px; white-space: nowrap;",
    "th": "border: 1px solid; border-width: 1px 0; text-align: left;",
}
styler.set_table_styles([dict(selector=sel, props=styles[sel]) for sel in styles])
styler.set_uuid("")
# hide redundant metrics (TPR = Recall, FPR = 1 - TNR, FNR = 1 - TPR)
styler.hide(["Recall", "FPR", "FNR"], axis=1)


# %% export model metrics as styled HTML table
# insert svelte {...props} forwarding to the table element
insert = """
<script>
  import { sortable } from 'svelte-zoo/actions'
</script>

<table use:sortable {...$$props}
"""
html_table = styler.to_html().replace("<table", insert)
with open(f"{FIGS}/metrics-table.svelte", "w") as file:
    file.write(html_table)


# %%
# hide_rows = list(set(df_metrics) - set(df_metrics.T.F1.nlargest(6).index))
# styler.hide(hide_rows)  # show only the best models by F1 score
dfi.export(styler, f"{ROOT}/tmp/figures/model-metrics.png", dpi=300)


# %% write model metrics to json for use by the website
df_stats["missing_preds"] = df_wbm[list(df_metrics)].isna().sum()
df_stats["missing_percent"] = [f"{x / len(df_wbm):.2%}" for x in df_stats.missing_preds]

df_stats.attrs["Total Run Time"] = df_stats[time_col].sum()

df_stats.round(2).to_json(f"{MODELS}/model-stats.json", orient="index")


# %% plot model run times as pie chart
fig = px.pie(df_stats, values=time_col, names=df_stats.index, hole=0.5).update_traces(
    textinfo="percent+label",
    textfont_size=14,
    marker=dict(line=dict(color="#000000", width=2)),
    hoverinfo="label+percent+name",
    texttemplate="%{label}<br>%{percent:.1%}",
    hovertemplate="%{label} %{percent:.1%} (%{value:.1f} h)",
    rotation=90,
    showlegend=False,
)
fig.add_annotation(
    # add title in the middle saying "Total CPU+GPU time used"
    text=f"Total CPU+GPU<br>time used:<br>{df_stats[time_col].sum():.1f} h",
    font=dict(size=18),
    x=0.5,
    y=0.5,
    showarrow=False,
)
fig.update_layout(margin=dict(l=0, r=0, t=0, b=0))


# %%
save_fig(fig, f"{FIGS}/model-run-times-pie.svelte")
