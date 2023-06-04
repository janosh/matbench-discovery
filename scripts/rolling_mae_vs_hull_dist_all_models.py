"""Plot rolling MAE as a function of hull distance for all models."""


# %%
from typing import Final

from pymatviz.utils import save_fig

from matbench_discovery import FIGS, PDF_FIGS
from matbench_discovery.plots import rolling_mae_vs_hull_dist
from matbench_discovery.preds import (
    df_each_pred,
    df_metrics,
    df_preds,
    each_true_col,
    models,
)

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

df_err, df_std = None, None  # variables to cache rolling MAE and std


# %%
# sort df columns by MAE (so that the legend is sorted too)
backend: Final = "plotly"

fig, df_err, df_std = rolling_mae_vs_hull_dist(
    e_above_hull_true=df_preds[each_true_col],
    e_above_hull_errors=df_each_pred[models],
    backend=backend,
    with_sem=False,
    df_rolling_err=df_err,
    df_err_std=df_std,
    show_dummy_mae=False,
)

if backend == "matplotlib":
    # increase line width in legend
    legend = fig.legend(frameon=False, loc="lower right")
    fig.figure.set_size_inches(10, 9)
    for handle in legend.get_lines():
        handle._linewidth *= 6
    for line in fig.lines:
        line._linewidth *= 2
else:
    for trace in fig.data:
        model = trace.name.split(" MAE=")[0]
        if model in df_metrics.T.sort_values("MAE").index[6:]:
            trace.visible = "legendonly"  # initially show only top models

    # increase line width
    fig.update_traces(line=dict(width=3))
    fig.layout.legend.update(bgcolor="rgba(0,0,0,0)")
    # increase legend handle size and reverse order
    fig.layout.margin.update(l=5, r=5, t=5, b=55)
    fig.show()


# %%
img_name = "rolling-mae-vs-hull-dist-models"
save_fig(fig, f"{FIGS}/{img_name}.svelte")
save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf", width=520, height=350)
