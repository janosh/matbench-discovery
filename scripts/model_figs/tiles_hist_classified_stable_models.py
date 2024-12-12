"""Histogram of the energy difference (either according to DFT ground truth [default] or
model predicted energy) to the convex hull for materials in the WBM data set. The
histogram stacks true/false positives/negatives with different colors.
"""

# %%
import math
from typing import Final

import pymatviz as pmv
from pymatviz.enums import Key

from matbench_discovery import PDF_FIGS, SITE_FIGS
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.models import MODEL_METADATA, model_is_compliant
from matbench_discovery.plots import hist_classified_stable_vs_hull_dist
from matbench_discovery.preds.discovery import (
    df_metrics,
    df_metrics_uniq_protos,
    df_preds,
    models,
)

__author__ = "Janosh Riebesell"
__date__ = "2022-12-01"


test_subset = globals().get("test_subset", TestSubset.uniq_protos)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(MbdKey.uniq_proto)
    df_metrics = df_metrics_uniq_protos

show_non_compliant = globals().get("show_non_compliant", False)
models_to_plot = [
    model
    for model in models
    if show_non_compliant or model_is_compliant(MODEL_METADATA[model])
]


# %%
n_cols = 3
use_full_rows = globals().get("use_full_rows", True)
if use_full_rows:
    # drop last models that don't fit in last row
    n_rows = len(models_to_plot) // n_cols
    models_to_plot = models_to_plot[: n_rows * n_cols]
else:
    n_rows = math.ceil(len(models) / n_cols)

hover_cols = (df_preds.index.name, MbdKey.e_form_dft, MbdKey.each_true, Key.formula)
facet_col = "Model"

# 'true' or 'pred': whether to put DFT or model-predicted hull distances on the x-axis
which_energy: Final = "pred"
kwds = dict(
    facet_col=facet_col,
    facet_col_wrap=n_cols,
    category_orders={facet_col: models_to_plot},
    facet_col_spacing=0.04,
    facet_row_spacing=0.04,
)

df_melt = df_preds.melt(
    id_vars=hover_cols,
    value_vars=models_to_plot,
    var_name=facet_col,
    value_name=Key.e_form_pred,
)

df_melt[Key.each_pred] = (
    df_melt[MbdKey.each_true] + df_melt[Key.e_form_pred] - df_melt[MbdKey.e_form_dft]
)

fig = hist_classified_stable_vs_hull_dist(
    df=df_melt,
    each_true_col=MbdKey.each_true,
    each_pred_col=Key.each_pred,
    which_energy=which_energy,
    rolling_acc=None,
    stability_threshold=None,
    **kwds,  # type: ignore[arg-type]
)

# TODO add line showing the true hull distance histogram on each subplot
show_metrics = False
for anno in fig.layout.annotations:
    model_name = anno.text = anno.text.split("=", 1).pop()
    if model_name not in models_to_plot or not show_metrics:
        continue
    F1, FPR, FNR, DAF = (df_metrics[model_name][x] for x in "F1 FPR FNR DAF".split())
    anno.text = f"{model_name} · {F1=:.2f} · {FPR=:.2f} · {FNR=:.2f} · {DAF=:.2f}"

# set the figure size based on the number of rows and columns
fig.layout.height = 230 * n_rows
fig.layout.width = 280 * n_cols

# set the shared y and x axis ranges
fig.update_yaxes(range=[0, 9_000], title_text=None, matches=None)
fig.update_xaxes(range=[-0.4, 0.4], title_text=None, matches=None)

axis_titles = dict(xref="paper", yref="paper", showarrow=False, font_size=16)
fig.add_annotation(  # x-axis title
    x=0.5,
    y=0,
    yshift=-50,
    text=MbdKey.each_true.label,
    borderpad=5,
    **axis_titles,
)
fig.add_annotation(  # y-axis title
    x=0,
    xshift=-70,
    y=0.5,
    text=Key.each_pred.label,
    textangle=-90,
    borderpad=5,
    **axis_titles,
)

# place the legend above the subplots
fig.layout.legend.update(
    y=1.08, xanchor="center", x=0.5, bgcolor="rgba(0,0,0,0)", orientation="h"
)

# standardize the margins and template
portrait = n_rows > n_cols
fig.layout.margin.update(l=60, r=10, t=0 if portrait else 10, b=60 if portrait else 10)
fig.layout.template = "pymatviz_white"

# for trace in fig.data:
#     # no need to store all 250k x values in plot, leads to 1.7 MB file,
#     # subsample every 10th point is enough to see the distribution
#     trace.x = trace.x[::10]

# increase height of figure
fig.show()


# %%
img_suffix = "" if show_non_compliant else "-only-compliant"
img_name = f"hist-clf-{which_energy}-hull-dist-models-{n_rows}x{n_cols}{img_suffix}"
pmv.save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf")
