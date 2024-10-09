"""Plot rolling MAE as a function of hull distance for a single model."""

# %%
import pymatviz as pmv
from pymatviz.enums import Key
from pymatviz.utils import MATPLOTLIB, PLOTLY

from matbench_discovery import PDF_FIGS, SITE_FIGS
from matbench_discovery.data import Model
from matbench_discovery.enums import MbdKey
from matbench_discovery.plots import rolling_mae_vs_hull_dist
from matbench_discovery.preds import df_each_pred, df_metrics, df_wbm

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"


# %%
model = Model.chgnet.label

fig, df_err, df_std = rolling_mae_vs_hull_dist(
    e_above_hull_true=df_wbm[MbdKey.each_true],
    e_above_hull_preds={model: df_each_pred[model]},
    backend=(backend := PLOTLY),
)

MAE, DAF, F1 = df_metrics[model][["MAE", Key.daf, "F1"]]
title = f"{model} · {MAE=:.2f} · {DAF=:.2f} · {F1=:.2f}"
if backend == MATPLOTLIB:
    fig = fig.figure
    fig.set_size_inches(6, 5)
    fig.legend(loc="lower right", frameon=False)
    fig.set(title=title)
    for line in fig.lines:
        line._linewidth *= 2  # noqa: SLF001
elif backend == PLOTLY:
    fig.update_layout(title=dict(text=title, x=0.5))
    fig.show()


# %%
img_name = f"rolling-mae-vs-hull-dist-{model}"
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf")
pmv.save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
