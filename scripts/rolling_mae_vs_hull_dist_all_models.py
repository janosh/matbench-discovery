# %%
from pymatviz.utils import save_fig

from matbench_discovery import FIGS, today
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.plots import Backend, rolling_mae_vs_hull_dist

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"


# %%
models = sorted(
    "Wrenformer, CGCNN, Voronoi RF, MEGNet, M3GNet, BOWSR MEGNet".split(", ")
)
e_form_col = "e_form_per_atom_mp2020_corrected"
e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"

df_wbm = load_df_wbm_with_preds(models=models).round(3)


# %%
backend: Backend = "plotly"

for model in models:
    model_error = df_wbm[e_form_col] - df_wbm[model]
    MAE = (df_wbm[e_above_hull_col] - model_error).abs().mean()
    df_wbm[f"{model} {MAE=:.2f}"] = df_wbm[e_form_col] - df_wbm[model]

fig, df_err, df_std = rolling_mae_vs_hull_dist(
    e_above_hull_true=df_wbm[e_above_hull_col],
    e_above_hull_errors=df_wbm.filter(like=" MAE="),
    backend=backend,
    # template="plotly_white",
)


if backend == "matplotlib":
    # increase line width in legend
    legend = fig.legend(frameon=False, loc="lower right")
    fig.figure.set_size_inches(10, 9)
    for handle in legend.get_lines():
        handle._linewidth *= 6
    for line in fig.lines:
        line._linewidth *= 3
else:
    fig.show()


# %%
img_path = f"{FIGS}/{today}-rolling-mae-vs-hull-dist-compare-models"
save_fig(fig, f"{img_path}.pdf")
