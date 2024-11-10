"""Plot ROC and PR (precision-recall) curves for each model."""

# %%
import math

import pandas as pd
import pymatviz as pmv
from pymatviz.enums import Key
from pymatviz.utils import PLOTLY
from sklearn.metrics import auc, roc_curve
from tqdm import tqdm

from matbench_discovery import PDF_FIGS, SITE_FIGS, STABILITY_THRESHOLD
from matbench_discovery import plots as plots
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.models import MODEL_METADATA, model_is_compliant
from matbench_discovery.preds.discovery import (
    df_each_pred,
    df_preds,
    model_styles,
    models,
)

__author__ = "Janosh Riebesell"
__date__ = "2023-01-30"


line = dict(dash="dash", width=0.5)
facet_col = "Model"
color_col = "Stability Threshold"
n_cols = 3
n_rows = math.ceil(len(models) / n_cols)


test_subset = globals().get("test_subset", TestSubset.uniq_protos)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(Key.uniq_proto)
    df_each_pred = df_each_pred.loc[df_preds.index]

show_non_compliant = globals().get("show_non_compliant", False)
models_to_plot = [
    model
    for model in models
    if show_non_compliant or model_is_compliant(MODEL_METADATA[model])
]
df_each_pred = df_each_pred[models_to_plot]


# %%
df_roc = pd.DataFrame()

for model in (pbar := tqdm(models_to_plot, desc="Calculating ROC curves")):
    pbar.set_postfix_str(model)

    na_mask = df_preds[MbdKey.each_true].isna() | df_each_pred[model].isna()
    y_true = (df_preds[~na_mask][MbdKey.each_true] <= STABILITY_THRESHOLD).astype(int)
    y_pred = df_each_pred[model][~na_mask]
    fpr, tpr, thresholds = roc_curve(y_true, y_pred, pos_label=0)
    AUC = auc(fpr, tpr)
    title = f"{model} · {AUC=:.2f}"
    thresholds = [f"{thresh:.3} eV/atom" for thresh in thresholds]
    df_tmp = pd.DataFrame(
        {"FPR": fpr, "TPR": tpr, color_col: thresholds, "AUC": AUC, facet_col: title}
    ).round(3)

    df_roc = pd.concat([df_roc, df_tmp])


# %%
facet_plot = False
kwds = dict(
    height=150 * len(df_roc[facet_col].unique()),
    color=color_col,
    facet_col=facet_col,
    range_color=(-0.5, 0.5),
    facet_col_spacing=0.03,
    facet_row_spacing=0.1,
)

plot_fn = getattr(
    df_roc.iloc[:: len(df_roc) // 500 or 1]
    .sort_values(["AUC", "FPR"], ascending=False)
    .plot,
    "scatter" if facet_plot else "line",
)

fig = plot_fn(
    x="FPR",
    y="TPR",
    facet_col_wrap=n_cols,
    backend=PLOTLY,
    range_x=(-0.01, 1),
    range_y=(0, 1.02),
    hover_name=facet_col,
    hover_data={facet_col: False, color_col: True},
    **(kwds if facet_plot else dict(color=facet_col, markers=True)),
)

for anno in fig.layout.annotations:
    anno.text = anno.text.split("=", 1)[1]  # remove Model= from subplot titles


for trace in fig.data:
    if styles := model_styles.get(trace.name.split(" · ")[0]):
        ls, marker, color = styles
        trace.line = dict(color=color, dash=ls, width=2)
        trace.marker = dict(color=color, symbol=marker, size=4)


if not facet_plot:
    fig.layout.legend.update(
        x=1,
        y=0,
        xanchor="right",
        title=None,
        bgcolor="rgba(0,0,0,0)",  # Transparent background
    )

fig.layout.coloraxis.colorbar.update(thickness=14, title_side="right")
if n_cols == 2:
    fig.layout.coloraxis.colorbar.update(
        x=1, y=1, xanchor="right", yanchor="top", lenmode="pixels", len=210
    )

fig.add_shape(type="line", x0=0, y0=0, x1=1, y1=1, line=line, row="all", col="all")
fig.add_annotation(
    text="No skill", x=0.5, y=0.5, showarrow=False, yshift=10, textangle=-45
)
# allow scrolling and zooming each subplot individually
fig.update_xaxes(matches=None)
fig.layout.margin.update(l=0, r=0, b=0, t=20, pad=0)
fig.update_yaxes(matches=None)
fig.show()
img_name = (
    f"roc-models{f'-{n_rows}x{n_cols}' if facet_plot else ''}"
    f"{'-only-compliant' if not show_non_compliant else ''}"
)


# %%
pmv.save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf", width=500, height=500)
