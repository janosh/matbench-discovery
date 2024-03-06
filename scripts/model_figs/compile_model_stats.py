"""Compile metrics and total run times for all models and export them to JSON, a
pandas-styled HTML table and a plotly figure.
"""

# %%
from __future__ import annotations

import re
from typing import Any

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import requests
import wandb
import wandb.apis.public
from pymatviz.io import save_fig
from tqdm import tqdm

from matbench_discovery import PDF_FIGS, SITE_FIGS, SITE_LIB, WANDB_PATH
from matbench_discovery.preds import (
    df_metrics,
    df_metrics_10k,
    df_metrics_uniq_protos,
    df_preds,
    model_styles,
)
from matbench_discovery.preds import models as all_models

__author__ = "Janosh Riebesell"
__date__ = "2022-11-28"


# %% wandb API queries for relevant model runs to compute train and test times from
train_run_filters: dict[str, tuple[int, str, str, str]] = {
    # model: (n_runs, display_name, created_gt, created_lt)
    "CGCNN": (10, "train-cgcnn-ensemble", "2022-11-21", "2022-11-23"),
    "Voronoi RF": (48, "voronoi-features-mp", "2022-11-24", "2022-11-26"),
    "Wrenformer": (10, "train-wrenformer-ensemble", "2022-11-14", "2022-11-16"),
}
test_run_filters: dict[str, tuple[int, str, str, str]] = {
    # model: (n_runs, display_name, created_gt, created_lt)
    "BOWSR": (476, "bowsr-megnet", "2023-01-20", "2023-01-22"),
    "CHGNet": (100, "chgnet-wbm-IS2RE-", "2023-03-05", "2023-03-07"),  # v0.2.0
    # "CHGNet": (100, "chgnet-wbm-IS2RE-", "2023-10-22", "2023-10-25"),  # v0.3.0
    "CGCNN": (1, "test-cgcnn-wbm-IS2RE", "2022-12-03", "2022-12-05"),
    "MACE": (100, "mace-wbm-IS2RE-FIRE", "2023-07-22", "2023-07-24"),
    "M3GNet": (99, "m3gnet-wbm-IS2RE", "2022-10-31", "2022-11-01"),
    "MEGNet": (1, "megnet-wbm-IS2RE", "2022-11-17", "2022-11-19"),
    "Voronoi RF": (20, "voronoi-features-wbm", "2022-11-15", "2022-11-19"),
    "Wrenformer": (1, "test-wrenformer", "2022-11-23", "2022-11-24"),
}

train_stats: dict[str, dict[str, Any]] = {}
test_stats: dict[str, dict[str, Any]] = {}


# %% calculate total model run times from wandb logs
# NOTE these model run times are pretty meaningless since some models were run on GPU
# (Wrenformer and CGCNN), others on CPU. Also BOWSR, M3GNet and MEGNet weren't
# trained from scratch. Their run times only indicate the time needed to predict the
# test set.

time_col = "Run Time (h)"
for label, stats, raw_filters in (
    ("train", train_stats, train_run_filters),
    ("test", test_stats, test_run_filters),
):
    for model in (pbar := tqdm(raw_filters, desc=f"Get WandB {label} runs")):
        n_runs, name, created_gt, created_lt = raw_filters[model]

        date_pat = r"\d{4}-\d{2}-\d{2}"
        assert re.match(date_pat, created_gt), f"{created_gt=} must be yyyy-mm-dd"
        assert re.match(date_pat, created_lt), f"{created_lt=} must be yyyy-mm-dd"
        assert created_gt < created_lt, f"{created_gt=} must be before {created_lt=}"

        if n_runs == 0 or model in stats:
            continue  # skip models with no runs or already processed

        pbar.set_postfix_str(model)
        filters = {  # construct wandb API query
            "display_name": {"$regex": name},
            "created_at": {"$gt": created_gt, "$lt": created_lt},
        }
        runs = wandb.Api().runs(WANDB_PATH, filters=filters)

        assert (
            len(runs) == n_runs
        ), f"found {len(runs)} {label} runs for {name!r}, expected {n_runs}"

        run_times = [run.summary.get("_wandb", {}).get("runtime", 0) for run in runs]

        run_time_total = sum(run_times)
        # NOTE we assume all jobs have the same metadata here
        metadata = requests.get(
            runs[0].file("wandb-metadata.json").url, timeout=10
        ).json()

        n_gpu, n_cpu = metadata.get("gpu_count", 0), metadata.get("cpu_count", 0)
        stats[model] = {
            time_col: run_time_total / 3600,
            "GPU": n_gpu,
            "CPU": n_cpu,
            "Slurm Jobs": n_runs,
        }

if False:  # commented out due to not incl. M3GNet + MEGNet in main analysis
    # add M3GNet and MEGNet run times for M3GNet + MEGNet
    test_stats["M3GNet + MEGNet"] = test_stats["M3GNet"].copy()
    test_stats["M3GNet + MEGNet"][time_col] = (
        test_stats["MEGNet"][time_col] + test_stats["M3GNet"][time_col]
    )
test_stats["CGCNN+P"] = {}


# %%
stats_dict = {
    "": df_metrics,
    "-10k": df_metrics_10k,
    "-uniq-protos": df_metrics_uniq_protos,
}
for label, df_tmp in stats_dict.items():
    df_train_stats = pd.DataFrame(train_stats).add_prefix("Train ", axis="index")
    df_test_stats = pd.DataFrame(test_stats).add_prefix("Test ", axis="index")
    df_tmp = pd.concat([df_tmp, df_train_stats, df_test_stats]).T

    df_tmp[time_col] = df_tmp.filter(like=time_col).sum(axis="columns")

    # write model metrics to json for website use
    df_tmp["missing_preds"] = df_preds[all_models].isna().sum()
    df_tmp["missing_percent"] = [
        f"{x / len(df_preds):.2%}" for x in df_tmp.missing_preds
    ]

    df_tmp.attrs["All Models Run Time"] = df_tmp[time_col].sum()

    # write stats for different data subsets to JSON
    df_tmp.round(2).to_json(f"{SITE_LIB}/model-stats{label}.json", orient="index")
    stats_dict[label] = df_tmp

df_stats = stats_dict[""]


# %%
df_time = (
    df_stats.sort_index()
    .filter(like=time_col)
    .round(1)
    # maybe remove BOWSR since it used so much more compute time than the other models
    # that it makes the plot unreadable
    # .drop(index="BOWSR")
    .reset_index(names=(model_col := "Model"))
)


# %% plot model run times as pie chart
fig = px.pie(
    df_time,
    values=time_col,
    names=model_col,
    hole=0.4,
    width=600,
    height=600,
).update_traces(
    textinfo="percent+label",
    textfont_size=14,
    marker=dict(line=dict(color="black", width=2)),
    hoverinfo="label+percent+name",
    texttemplate="%{label}<br>%{percent:.1%}",
    hovertemplate="%{label} %{percent:.1%} (%{value:.1f} h)",
    rotation=90,
    showlegend=False,
)
fig.layout.margin.update(l=0, r=0, t=0, b=0)
title = f"Total CPU+GPU<br>time used:<br>{df_stats[time_col].sum():.1f} h"
fig.add_annotation(text=title, font=dict(size=15), x=0.5, y=0.5, showarrow=False)
pie_path = f"{SITE_FIGS}/model-run-times-pie.svelte"
# save_fig(fig, pie_path)
fig.show()


# %% plot model run times as sunburst chart
run_type_col = "Run Type"
df_melt = (
    df_time.rename(columns=lambda col: col.replace(f" {time_col}", ""))
    .melt(id_vars=model_col, value_vars=["Train", "Test"], var_name=run_type_col)
    .rename(columns={"value": time_col})
)
# see https://stackoverflow.com/a/69799623 for adding a hole
fig = px.sunburst(
    df_melt[df_melt[time_col] > 0],
    path=[model_col, run_type_col],
    values=time_col,
    color=model_col,
)
fig.layout.margin.update(l=0, r=0, t=0, b=0)
fig.update_traces(marker=dict(line=dict(color="white", width=1)))


# %% plot model run times as bar chart
fig = df_melt.dropna().plot.bar(
    x=time_col,
    y=model_col,
    backend="plotly",
    # color=time_col,
    text_auto=".0f",
    text=time_col,
    color=model_col,
    color_discrete_sequence=[model_styles[model][2] for model in df_melt[model_col]],
)
# reduce bar width
fig.update_traces(width=0.8)

title = f"All models: {df_stats[time_col].sum():.0f} h"
fig.layout.legend.update(title=title, orientation="h", xanchor="center", x=0.4, y=1.2)
fig.layout.xaxis.title = ""
fig.layout.margin.update(l=0, r=0, t=0, b=0)
fig.show()


# %%
save_fig(fig, f"{SITE_FIGS}/model-run-times-bar.svelte")
pdf_fig = go.Figure(fig)
# replace legend with annotation in PDF
pdf_fig.layout.showlegend = False
pdf_fig.add_annotation(
    text=title,
    font=dict(size=15),
    x=0.99,
    y=0.99,
    showarrow=False,
    xref="paper",
    yref="paper",
)
save_fig(pdf_fig, f"{PDF_FIGS}/model-run-times-bar.pdf", height=300, width=800)
