# %%
from typing import Final

from pymatviz.utils import save_fig

from matbench_discovery import FIGS, ROOT
from matbench_discovery.plots import rolling_mae_vs_hull_dist
from matbench_discovery.preds import df_each_pred, df_metrics, df_wbm, each_true_col

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"


# %%
# sort df columns by MAE (so that the legend is sorted too)
backend: Final = "plotly"

fig, df_err, df_std = rolling_mae_vs_hull_dist(
    e_above_hull_true=df_wbm[each_true_col],
    e_above_hull_errors=df_each_pred,
    backend=backend,
    with_sem=False,
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
    # keep only every n-th point to reduce plot size for website
    for trace in fig.data:
        if trace.name and trace.name.startswith("MAE") and len(trace.x) < 100:
            continue  # skip the MAE < DFT error area traces
        trace.x = trace.x[::5]
        trace.y = trace.y[::5]

        # initially show only the top 5 models
        if trace.name.split(" MAE=")[0] not in df_metrics.T.MAE.nsmallest(5):
            trace.visible = "legendonly"
            # trace.visible = False

    # increase line width
    fig.update_traces(line=dict(width=3))

    # increase legend handle size and reverse order
    fig.layout.margin = dict(l=5, r=5, t=5, b=55)
    fig.layout.legend.update(itemsizing="constant", bgcolor="rgba(0,0,0,0)")
    fig.show()


# %%
img_name = "rolling-mae-vs-hull-dist-models"
save_fig(fig, f"{FIGS}/{img_name}.svelte")
# save_fig(fig, f"{STATIC}/{img_name}.webp", scale=3)
save_fig(fig, f"{ROOT}/tmp/figures/{img_name}.pdf")
