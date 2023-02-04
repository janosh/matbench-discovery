# %%
from pymatviz.utils import save_fig

from matbench_discovery import FIGS, STATIC, today
from matbench_discovery.plots import Backend, rolling_mae_vs_hull_dist
from matbench_discovery.preds import df_metrics, df_wbm, e_form_col, each_true_col

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"


# %%
# sort df columns by MAE (so that the legend is sorted too)
backend: Backend = "plotly"

for model, MAE in sorted(df_metrics.T.MAE.items(), key=lambda x: x[1]):
    df_wbm[f"{model} {MAE=:.2f}"] = df_wbm[e_form_col] - df_wbm[model]

fig, df_err, df_std = rolling_mae_vs_hull_dist(
    e_above_hull_true=df_wbm[each_true_col],
    e_above_hull_errors=df_wbm.filter(like=" MAE="),
    backend=backend,
    with_sem=False,
    template="plotly_white" if (mode := "light") == "light" else None,
    width=800,
    height=800,
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

    # increase line width
    fig.update_traces(line=dict(width=3))

    # increase legend handle size and reverse order
    fig.layout.legend.update(itemsizing="constant", traceorder="reversed")
    fig.show()


# %%
img_path = f"{today}-rolling-mae-vs-hull-dist-models-{mode}"
save_fig(fig, f"{FIGS}/{img_path}.svelte")
save_fig(fig, f"{STATIC}/{img_path}.webp", scale=3)
