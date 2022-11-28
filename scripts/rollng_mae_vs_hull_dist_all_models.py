# %%
from matbench_discovery import ROOT, today
from matbench_discovery.load_preds import load_df_wbm_with_preds
from matbench_discovery.plots import plt, rolling_mae_vs_hull_dist

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"


# %%
models = (
    "Wren, CGCNN IS2RE, CGCNN RS2RE, Voronoi RF, "
    "Wrenformer, MEGNet, M3GNet, BOWSR MEGNet"
).split(", ")

df_wbm = load_df_wbm_with_preds(models=models).round(3)


# %%
fig, ax = plt.subplots(1, figsize=(10, 9))

target_col = "e_form_per_atom_mp2020_corrected"
e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"

for model_name in sorted(models):

    # assert df_wbm[model_name].isna().sum() < 100

    rolling_mae_vs_hull_dist(
        e_above_hull_pred=df_wbm[model_name] - df_wbm[target_col],
        e_above_hull_true=df_wbm[e_above_hull_col],
        ax=ax,
        label=model_name,
    )

# increase line width in legend
legend = ax.legend(frameon=False, loc="lower right")
for line in legend.get_lines():
    line._linewidth *= 3


# %%
img_path = f"{ROOT}/figures/{today}-rolling-mae-vs-hull-dist-compare-models.pdf"
# fig.savefig(img_path)
