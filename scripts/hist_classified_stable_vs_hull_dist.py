"""Histogram of the energy difference (either according to DFT ground truth [default] or
model predicted energy) to the convex hull for materials in the WBM data set. The
histogram stacks true/false positives/negatives with different colors.

See fig. S1 in https://science.org/doi/10.1126/sciadv.abn4117.
"""


# %%
from typing import Final

from pymatviz.io import save_fig

from matbench_discovery import PDF_FIGS
from matbench_discovery.data import df_wbm
from matbench_discovery.plots import hist_classified_stable_vs_hull_dist
from matbench_discovery.preds import df_each_pred, each_true_col

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"


# %%
# model_name = "Wrenformer"
model_name = "CGCNN"
# model_name = "CGCNN+P"
which_energy: Final = "pred"
df_each_pred[each_true_col] = df_wbm[each_true_col]
backend: Final = "plotly"

fig = hist_classified_stable_vs_hull_dist(
    df_each_pred,
    each_true_col=each_true_col,
    each_pred_col=model_name,
    which_energy=which_energy,
    # stability_threshold=-0.05,
    rolling_acc=None,
    backend=backend,
)

if backend == "plotly":
    # fig.layout.title.update(text=model_name, x=0.5)
    fig.layout.margin.update(l=0, r=0, b=0, t=30)
    # fig.update_yaxes(range=[0, 12000])
    fig.show()


# %%
img_name = f"hist-clf-{which_energy}-hull-dist-{model_name}"
# save_fig(fig, f"{FIGS}/{img_name}.svelte")
save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf")
