# %%
from plotly.subplots import make_subplots
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

df_wbm = load_df_wbm_with_preds(models=models).round(3)


# %%
target_col = "e_form_per_atom_mp2020_corrected"
e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"
backend: Backend = "plotly"

rows, cols = len(models) // 3, 3
if backend == "plotly":
    fig = make_subplots(rows=rows, cols=cols)


for idx, model_name in enumerate(models):
    row, col = idx % rows + 1, idx // rows + 1

    # assert df_wbm[model_name].isna().sum() < 100
    preds = df_wbm[target_col] - df_wbm[model_name]
    MAE = (df_wbm[e_above_hull_col] - preds).abs().mean()

    ax = rolling_mae_vs_hull_dist(
        e_above_hull_true=df_wbm[e_above_hull_col],
        e_above_hull_error=preds,
        label=f"{model_name} · {MAE=:.2f}",
        backend=backend,
    )
    if backend == "plotly":
        fig.add_traces(ax.data, row=row, col=col)

if hasattr(ax, "legend"):
    # increase line width in legend
    legend = ax.legend(frameon=False, loc="lower right")
    ax.figure.set_size_inches(10, 9)
    for line in legend.get_lines():
        line._linewidth *= 3


fig.show()


# %%
img_path = f"{FIGS}/{today}-rolling-mae-vs-hull-dist-compare-models"
save_fig(fig, f"{img_path}.pdf")
