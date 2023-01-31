# %%
from pymatviz.utils import save_fig

from matbench_discovery import FIGS, today
from matbench_discovery.data import load_df_wbm_preds
from matbench_discovery.plots import Backend, rolling_mae_vs_hull_dist

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"


# %%
models = sorted(
    "Wrenformer, CGCNN, Voronoi Random Forest, MEGNet, M3GNet, BOWSR MEGNet".split(", ")
)
e_form_col = "e_form_per_atom_mp2020_corrected"
e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"

df_wbm = load_df_wbm_preds(models=models).round(3)


# %%
backend: Backend = "plotly"


MAEs = {}
for model in models:
    # e_form and e_above_hull MAE are the same so we compute former for simplicity and
    # use it in place of the latter
    MAE = (df_wbm[e_form_col] - df_wbm[model]).abs().mean()
    MAEs[MAE] = model

# sort df columns by MAE (so that the legend is sorted too)
for MAE, model in sorted(MAEs.items(), reverse=True):
    df_wbm[f"{model} {MAE=:.2f}"] = df_wbm[e_form_col] - df_wbm[model]

fig, df_err, df_std = rolling_mae_vs_hull_dist(
    e_above_hull_true=df_wbm[e_above_hull_col],
    e_above_hull_errors=df_wbm.filter(like=" MAE="),
    backend=backend,
    with_sem=False,
    # template="plotly_white",
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
    # keep every 10th point from each trace to reduce plot size for website
    for trace in fig.data:
        if trace.name and trace.name.startswith("MAE") and len(trace.x) < 100:
            continue  # skip the MAE < DFT error area traces
        trace.x = trace.x[::10]
        trace.y = trace.y[::10]

    # increase line width
    fig.update_traces(line=dict(width=3))

    # increase legend handle size and reverse order
    fig.update_layout(legend=dict(itemsizing="constant"), legend_traceorder="reversed")

    fig.show()


# %%
img_path = f"{FIGS}/{today}-rolling-mae-vs-hull-dist-models"
save_fig(fig, f"{img_path}.svelte")
