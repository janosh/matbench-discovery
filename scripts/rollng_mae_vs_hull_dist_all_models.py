# %%
from matbench_discovery import FIGS, today
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.plots import rolling_mae_vs_hull_dist

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"


# %%
models = (
    "Wren, CGCNN, CGCNN IS2RE, CGCNN RS2RE, Voronoi RF, "
    "Wrenformer, MEGNet, M3GNet, BOWSR MEGNet"
).split(", ")

df_wbm = load_df_wbm_with_preds(models=models).round(3)


# %%
target_col = "e_form_per_atom_mp2020_corrected"
e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"

for model_name in sorted(models):

    # assert df_wbm[model_name].isna().sum() < 100
    preds = df_wbm[target_col] - df_wbm[model_name]
    MAE = (df_wbm[e_above_hull_col] - preds).abs().mean()

    ax = rolling_mae_vs_hull_dist(
        e_above_hull_true=df_wbm[e_above_hull_col],
        e_above_hull_error=preds,
        label=f"{model_name} · {MAE=:.2f}",
        backend="matplotlib",
    )

# increase line width in legend
legend = ax.legend(frameon=False, loc="lower right")
ax.figure.set_size_inches(10, 9)
for line in legend.get_lines():
    line._linewidth *= 3


# %%
img_path = f"{FIGS}/{today}-rolling-mae-vs-hull-dist-compare-models"
# ax.figure.savefig(f"{img_path}.pdf")
