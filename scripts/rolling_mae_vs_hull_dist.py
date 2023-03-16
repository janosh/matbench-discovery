"""Plot rolling MAE as a function of hull distance for a single model."""


# %%
from matbench_discovery import FIGS, today
from matbench_discovery.plots import rolling_mae_vs_hull_dist
from matbench_discovery.preds import df_metrics, df_preds, e_form_col, each_true_col

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"


# %%
model = "Wrenformer"
model = "MEGNet"
model = "CHGNet"

ax, df_err, df_std = rolling_mae_vs_hull_dist(
    e_above_hull_true=df_preds[each_true_col],
    e_above_hull_errors={model: df_preds[e_form_col] - df_preds[model]},
    # label=model,
    backend=(backend := "plotly"),
    # template="plotly_white",
)

MAE, DAF, F1 = df_metrics[model][["MAE", "DAF", "F1"]]
title = f"{today} {model} · {MAE=:.2f} · {DAF=:.2f} · {F1=:.2f}"
if backend == "matplotlib":
    fig = ax.figure
    fig.set_size_inches(6, 5)
    ax.legend(loc="lower right", frameon=False)
    ax.set(title=title)
    for line in ax.lines:
        line._linewidth *= 2
elif backend == "plotly":
    ax.update_layout(title=dict(text=title, x=0.5))
    ax.show()

img_path = f"{FIGS}/rolling-mae-vs-hull-dist.pdf"
# fig.savefig(img_path)
