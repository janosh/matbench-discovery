"""Compile metrics and total run times for all models and export them to JSON, a
pandas-styled HTML table and a plotly figure.
"""


# %%
from __future__ import annotations

import re
from typing import Any

import dataframe_image as dfi
import pandas as pd
import requests
import wandb
import wandb.apis.public
from pymatviz.utils import save_fig
from tqdm import tqdm

from matbench_discovery import FIGS, ROOT, WANDB_PATH
from matbench_discovery.plots import px
from matbench_discovery.preds import df_metrics, df_preds

__author__ = "Janosh Riebesell"
__date__ = "2022-11-28"


# %%
train_run_filters: dict[str, tuple[int, str, str, str]] = {
    # model: (n_runs, display_name, created_gt, created_lt)
    "CGCNN": (10, "train-cgcnn-ensemble", "2022-11-21", "2022-11-23"),
    "Voronoi RF": (48, "voronoi-features-mp", "2022-11-24", "2022-11-26"),
    "Wrenformer": (10, "train-wrenformer-ensemble", "2022-11-14", "2022-11-16"),
}
test_run_filters: dict[str, tuple[int, str, str, str]] = {
    # model: (n_runs, display_name, created_gt, created_lt)
    "BOWSR + MEGNet": (476, "bowsr-megnet", "2023-01-20", "2023-01-22"),
    "CHGNet": (100, "chgnet-wbm-IS2RE-", "2023-03-05", "2023-03-07"),
    "CGCNN": (1, "test-cgcnn-wbm-IS2RE", "2022-12-03", "2022-12-05"),
    "M3GNet": (99, "m3gnet-wbm-IS2RE", "2022-10-31", "2022-11-01"),
    "MEGNet": (1, "megnet-wbm-IS2RE", "2022-11-17", "2022-11-19"),
    "Voronoi RF": (20, "voronoi-features-wbm", "2022-11-15", "2022-11-19"),
    "Wrenformer": (1, "test-wrenformer", "2022-11-23", "2022-11-24"),
}

train_stats: dict[str, dict[str, Any]] = {}
test_stats: dict[str, dict[str, Any]] = {}


# %% calculate total model run times from wandb logs
# NOTE these model run times are pretty meaningless since some models were run on GPU
# (Wrenformer and CGCNN), others on CPU. Also BOWSR + MEGNet, M3GNet and MEGNet weren't
# trained from scratch. Their run times only indicate the time needed to predict the
# test set.

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
            continue
        pbar.set_postfix_str(model)
        filters = {
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
        metadata = requests.get(runs[0].file("wandb-metadata.json").url).json()

        n_gpu, n_cpu = metadata.get("gpu_count", 0), metadata.get("cpu_count", 0)
        stats[model] = {
            (time_col := "Run Time (h)"): run_time_total / 3600,
            "GPU": n_gpu,
            "CPU": n_cpu,
            "Slurm Jobs": n_runs,
        }

# test_stats["M3GNet + MEGNet"] = test_stats["M3GNet"].copy()
# test_stats["M3GNet + MEGNet"][time_col] = (
#     test_stats["MEGNet"][time_col] + test_stats["M3GNet"][time_col]
# )
test_stats["CGCNN+P"] = {}


df_stats = pd.concat(
    [
        df_metrics,
        pd.DataFrame(train_stats).add_prefix("Train ", axis="index"),
        pd.DataFrame(test_stats).add_prefix("Test ", axis="index"),
    ],
).T
df_stats[time_col] = df_stats.filter(like=time_col).sum(axis="columns")


# %%
time_cols = list(df_stats.filter(like=time_col))
for col in time_cols:  # uncomment to include run times in metrics table
    df_metrics.loc[col] = df_stats[col]
higher_is_better = {"DAF", "R²", "Precision", "F1", "Accuracy", "TPR", "TNR"}
lower_is_better = {"MAE", "RMSE", "FNR", "FPR", *time_cols}
idx_set = set(df_metrics.index)
styler = (
    df_metrics.T.rename(columns={"R2": "R²"})
    # append arrow up/down to table headers to indicate higher/lower metric is better
    # .rename(columns=lambda x: x + " ↑" if x in higher_is_better else x + " ↓")
    .style.format(precision=2)
    # reverse color map if lower=better
    .background_gradient(cmap="viridis_r", subset=list(lower_is_better & idx_set))
    # .background_gradient(
    #     cmap="viridis_r",
    #     subset=[time_col],
    #     gmap=np.log10(df_stats[time_col].to_numpy()),  # for log scaled color map
    # )
    .background_gradient(cmap="viridis", subset=list(higher_is_better & idx_set))
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
png_metrics = f"{ROOT}/tmp/figures/model-metrics.png"
dfi.export(styler, png_metrics, dpi=300)
print(f"{png_metrics=}")


# %% write model metrics to json for use by the website
df_stats["missing_preds"] = df_preds[list(df_metrics)].isna().sum()
df_stats["missing_percent"] = [
    f"{x / len(df_preds):.2%}" for x in df_stats.missing_preds
]

df_stats.attrs["All Models Run Time"] = df_stats[time_col].sum()
print(f"{df_stats[time_col].sum()=:.0f} hours")

# df_stats.round(2).to_json(f"{MODELS}/model-stats.json", orient="index")


# %% plot model run times as pie chart
df_time = (
    df_stats.sort_index()
    .filter(like=time_col)
    .round(1)
    # remove BOWSR since it used so much more compute time than the other models that
    # it makes the plot unreadable TODO decide if we want to keep BOWSR
    # .drop(index="BOWSR + MEGNet")
    .reset_index(names=(model_col := "Model"))
)
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
    marker=dict(line=dict(color="#000000", width=2)),
    hoverinfo="label+percent+name",
    texttemplate="%{label}<br>%{percent:.1%}",
    hovertemplate="%{label} %{percent:.1%} (%{value:.1f} h)",
    rotation=90,
    showlegend=False,
)
fig.layout.margin.update(l=0, r=0, t=0, b=0)
title = f"Total CPU+GPU<br>time used:<br>{df_stats[time_col].sum():.1f} h"
fig.add_annotation(text=title, font=dict(size=15), x=0.5, y=0.5, showarrow=False)
pie_path = f"{FIGS}/model-run-times-pie.svelte"
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
    y=time_col,
    x=model_col,
    backend="plotly",
    # color=time_col,
    text_auto=".0f",
    text=time_col,
    color=model_col,
)
legend_title = "Click to toggle models"
fig.layout.legend.update(
    x=0.98, y=0.98, xanchor="right", yanchor="top", title=legend_title
)
fig.layout.xaxis.title = ""
fig.layout.margin.update(l=0, r=0, t=0, b=0)
bar_path = f"{FIGS}/model-run-times-bar.svelte"
save_fig(fig, bar_path)
fig.show()
