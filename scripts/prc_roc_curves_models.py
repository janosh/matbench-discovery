"""Histogram of the energy difference (either according to DFT ground truth [default] or
model predicted energy) to the convex hull for materials in the WBM data set. The
histogram stacks true/false positives/negatives with different colors.
"""


# %%
import pandas as pd
from pymatviz.utils import save_fig
from sklearn.metrics import auc, precision_recall_curve, roc_curve
from tqdm import tqdm

from matbench_discovery import FIGS
from matbench_discovery.plots import pio
from matbench_discovery.preds import df_each_pred, df_wbm, each_true_col

__author__ = "Janosh Riebesell"
__date__ = "2023-01-30"


pio.templates.default
line = dict(dash="dash", width=0.5)

facet_col = "Model"
color_col = "Stability Threshold"


# %%
df_roc = pd.DataFrame()

for model in (pbar := tqdm(list(df_each_pred), desc="Calculating ROC curves")):
    pbar.set_postfix_str(model)
    na_mask = df_wbm[each_true_col].isna() | df_each_pred[model].isna()
    y_true = (df_wbm[~na_mask][each_true_col] <= 0).astype(int)
    y_pred = df_each_pred[model][~na_mask]
    fpr, tpr, thresholds = roc_curve(y_true, y_pred, pos_label=0)
    AUC = auc(fpr, tpr)
    title = f"{model} · {AUC=:.2f}"
    df_tmp = pd.DataFrame(
        {"FPR": fpr, "TPR": tpr, color_col: thresholds, "AUC": AUC, facet_col: title}
    ).round(3)

    df_roc = pd.concat([df_roc, df_tmp])


# %%
fig = (
    df_roc.iloc[:: len(df_roc) // 500 or 1]
    .sort_values(["AUC", "FPR"], ascending=False)
    .plot.scatter(
        x="FPR",
        y="TPR",
        facet_col=facet_col,
        facet_col_wrap=2,
        backend="plotly",
        height=150 * len(df_roc[facet_col].unique()),
        color=color_col,
        range_x=(0, 1),
        range_y=(0, 1),
        range_color=(-0.5, 0.5),
        hover_name=facet_col,
        hover_data={facet_col: False},
    )
)

for anno in fig.layout.annotations:
    anno.text = anno.text.split("=", 1)[1]  # remove Model= from subplot titles

fig.layout.coloraxis.colorbar.update(
    x=1,
    y=1,
    xanchor="right",
    yanchor="top",
    thickness=14,
    lenmode="pixels",
    len=210,
    title_side="right",
)
fig.add_shape(type="line", x0=0, y0=0, x1=1, y1=1, line=line, row="all", col="all")
fig.add_annotation(text="No skill", x=0.5, y=0.5, showarrow=False, yshift=-10)
# allow scrolling and zooming each subplot individually
fig.update_xaxes(matches=None)
fig.update_yaxes(matches=None)
fig.show()


# %%
save_fig(fig, f"{FIGS}/roc-models.svelte")


# %%
df_prc = pd.DataFrame()

for model in (pbar := tqdm(list(df_each_pred), desc="Calculating ROC curves")):
    pbar.set_postfix_str(model)
    na_mask = df_wbm[each_true_col].isna() | df_each_pred[model].isna()
    y_true = (df_wbm[~na_mask][each_true_col] <= 0).astype(int)
    y_pred = df_each_pred[model][~na_mask]
    prec, recall, thresholds = precision_recall_curve(y_true, y_pred, pos_label=0)
    df_tmp = pd.DataFrame(
        {
            "Precision": prec[:-1],
            "Recall": recall[:-1],
            color_col: thresholds,
            facet_col: model,
        }
    ).round(3)

    df_prc = pd.concat([df_prc, df_tmp])


# %%
fig = df_prc.iloc[:: len(df_roc) // 500 or 1].plot.scatter(
    x="Recall",
    y="Precision",
    facet_col=facet_col,
    facet_col_wrap=2,
    backend="plotly",
    height=150 * len(df_roc[facet_col].unique()),
    color=color_col,
    range_x=(0, 1),
    range_y=(0.5, 1),
    range_color=(-0.5, 1),
    hover_name=facet_col,
    hover_data={facet_col: False},
)

for anno in fig.layout.annotations:
    anno.text = anno.text.split("=", 1)[1]  # remove Model= from subplot titles

fig.layout.coloraxis.colorbar.update(
    x=0.5, y=1.1, thickness=14, len=0.4, orientation="h"
)
fig.add_hline(y=0.5, line=line)
fig.add_annotation(
    text="No skill", x=0, y=0.5, showarrow=False, xanchor="left", xshift=10, yshift=10
)
# allow scrolling and zooming each subplot individually
fig.update_xaxes(matches=None)
fig.update_yaxes(matches=None)
fig.show()


# %%
save_fig(fig, f"{FIGS}/prc-models.svelte")
fig.update_yaxes(matches=None)
fig.show()
