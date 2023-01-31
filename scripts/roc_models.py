# %%
import numpy as np
import pandas as pd
from pymatviz.utils import save_fig

from matbench_discovery import FIGS, today
from matbench_discovery.data import load_df_wbm_preds
from matbench_discovery.energy import stable_metrics
from matbench_discovery.plots import pio

__author__ = "Janosh Riebesell"
__date__ = "2023-01-30"

"""
Histogram of the energy difference (either according to DFT ground truth [default] or
model predicted energy) to the convex hull for materials in the WBM data set. The
histogram stacks true/false positives/negatives with different colors.
"""

pio.templates.default


# %%
models = sorted(
    "CGCNN, Voronoi Random Forest, Wrenformer, MEGNet, M3GNet, BOWSR MEGNet".split(", ")
)
df_wbm = load_df_wbm_preds(models).round(3)

e_form_col = "e_form_per_atom_mp2020_corrected"
each_true_col = "e_above_hull_mp2020_corrected_ppd_mp"
each_pred_col = "e_above_hull_pred"
facet_col = "Model"
color_col = "Stability Threshold"


# %%
df_roc = pd.DataFrame()

for model in models:
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
    x=1,
    y=1,
    xanchor="right",
    yanchor="top",
    thickness=14,
    len=0.27,
    title_side="right",
)
fig.show()


# %%
save_fig(fig, f"{FIGS}/{today}-roc-models.svelte")
