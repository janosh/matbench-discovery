"""Histogram of the energy difference (either according to DFT ground truth [default] or
model predicted energy) to the convex hull for materials in the WBM data set. The
histogram stacks true/false positives/negatives with different colors.

See fig. S1 in https://science.org/doi/10.1126/sciadv.abn4117.
"""


# %%
from typing import Final

from pymatviz.utils import save_fig

from matbench_discovery import FIGS
from matbench_discovery.data import df_wbm
from matbench_discovery.plots import hist_classified_stable_vs_hull_dist
from matbench_discovery.preds import df_each_pred, each_true_col

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"


# %%
model_name = "Wrenformer"
model_name = "CHGNet"
# model_name = "M3GNet"
# model_name = "Voronoi RF"
which_energy: Final = "true"
df_each_pred[each_true_col] = df_wbm[each_true_col]
backend: Final = "plotly"

fig = hist_classified_stable_vs_hull_dist(
    df_each_pred,
    each_true_col=each_true_col,
    each_pred_col=model_name,
    which_energy=which_energy,
    # stability_threshold=-0.05,
    # rolling_acc=None,
    backend=backend,
)

if backend == "plotly":
    fig.layout.title = model_name
    fig.show()


# %%
img_path = f"{FIGS}/hist-clf-{which_energy}-hull-dist-{model_name}"
# save_fig(fig, f"{img_path}.svelte")
save_fig(fig, f"{img_path}.webp")
