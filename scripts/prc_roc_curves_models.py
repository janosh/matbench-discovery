# %%
import numpy as np
import pandas as pd
from pymatviz.utils import save_fig
from tqdm import tqdm

from matbench_discovery import FIGS, today
from matbench_discovery.metrics import (
    df_wbm,
    e_form_col,
    each_pred_col,
    each_true_col,
    models,
    stable_metrics,
)
from matbench_discovery.plots import pio

__author__ = "Janosh Riebesell"
__date__ = "2023-01-30"

"""
Histogram of the energy difference (either according to DFT ground truth [default] or
model predicted energy) to the convex hull for materials in the WBM data set. The
histogram stacks true/false positives/negatives with different colors.
"""

pio.templates.default
line = dict(dash="dash", width=0.5)

facet_col = "Model"
color_col = "Stability Threshold"


# %%
df_roc = pd.DataFrame()

for model in (pbar := tqdm(models)):
    pbar.set_description(model)
    df_wbm[f"{model}_{each_pred_col}"] = df_wbm[each_true_col] + (
        df_wbm[model] - df_wbm[e_form_col]
    )
    for stab_treshold in np.arange(-0.4, 0.4, 0.01):
        metrics = stable_metrics(
            df_wbm[each_true_col], df_wbm[f"{model}_{each_pred_col}"], stab_treshold
        )
        df_tmp = pd.DataFrame(
            {facet_col: model, color_col: stab_treshold, **metrics}, index=[0]
        )
        df_roc = pd.concat([df_roc, df_tmp])


df_roc = df_roc.round(3)


# %%
fig = df_roc.plot.scatter(
    x="FPR",
    y="TPR",
    facet_col=facet_col,
    facet_col_wrap=2,
    backend="plotly",
    height=800,
    color=color_col,
    range_x=(0, 1),
    range_y=(0, 1),
)

for anno in fig.layout.annotations:
    anno.text = anno.text.split("=")[1]  # remove Model= from subplot titles

fig.layout.coloraxis.colorbar.update(
    x=1, y=1, xanchor="right", yanchor="top", thickness=14, len=0.27, title_side="right"
)
fig.add_shape(type="line", x0=0, y0=0, x1=1, y1=1, line=line, row="all", col="all")
fig.add_annotation(
    text="No skill", x=0.5, y=0.5, showarrow=False, yshift=-10, textangle=-30
)
# allow scrolling and zooming each subplot individually
fig.update_xaxes(matches=None)
fig.update_yaxes(matches=None)
fig.show()


# %%
save_fig(fig, f"{FIGS}/{today}-roc-models.svelte")


# %%
fig = df_roc.plot.scatter(
    x="Recall",
    y="Precision",
    facet_col=facet_col,
    facet_col_wrap=2,
    backend="plotly",
    height=800,
    color=color_col,
    range_x=(0, 1),
    range_y=(0, 1),
)

for anno in fig.layout.annotations:
    anno.text = anno.text.split("=")[1]  # remove Model= from subplot titles

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
save_fig(fig, f"{FIGS}/{today}-prc-models.svelte")
fig.update_yaxes(matches=None)
fig.show()
