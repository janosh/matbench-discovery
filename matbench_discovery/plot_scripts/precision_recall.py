# %%
from sklearn.metrics import f1_score

from matbench_discovery import ROOT, today
from matbench_discovery.plot_scripts import load_df_wbm_with_preds
from matbench_discovery.plots import StabilityCriterion, cumulative_clf_metric, plt

__author__ = "Rhys Goodall, Janosh Riebesell"


# %%
models = (
    "Wren, CGCNN IS2RE, CGCNN RS2RE, Voronoi RF, "
    "Wrenformer, MEGNet, M3GNet, BOWSR MEGNet"
).split(", ")

df_wbm = load_df_wbm_with_preds(models=models).round(3)

target_col = "e_form_per_atom_mp2020_corrected"
e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"


# %%
stability_crit: StabilityCriterion = "energy"
colors = "tab:blue tab:orange teal tab:pink black red turquoise tab:purple".split()


# %%
fig, (ax_prec, ax_recall) = plt.subplots(1, 2, figsize=(15, 7), sharey=True)

for model_name, color in zip(models, colors):

    e_above_hull_pred = df_wbm[model_name] - df_wbm[target_col]

    F1 = f1_score(df_wbm[e_above_hull_col] < 0, e_above_hull_pred < 0)

    e_above_hull_error = e_above_hull_pred + df_wbm[e_above_hull_col]
    cumulative_clf_metric(
        e_above_hull_error,
        df_wbm[e_above_hull_col],
        color=color,
        label=f"{model_name}\n{F1=:.3}",
        project_end_point="xy",
        stability_crit=stability_crit,
        ax=ax_prec,
        metric="precision",
    )

    cumulative_clf_metric(
        e_above_hull_error,
        df_wbm[e_above_hull_col],
        color=color,
        label=f"{model_name}\n{F1=:.3}",
        project_end_point="xy",
        stability_crit=stability_crit,
        ax=ax_recall,
        metric="recall",
    )


for ax in (ax_prec, ax_recall):
    ax.set(xlim=(0, None))


# x-ticks every 10k materials
# ax.set(xticks=range(0, int(ax.get_xlim()[1]), 10_000))

fig.suptitle(f"{today} {stability_crit=}")
xlabel_cumulative = "Materials predicted stable sorted by hull distance"
fig.text(0.5, -0.08, xlabel_cumulative, ha="center")


# %%
img_path = f"{ROOT}/figures/{today}-precision-recall-curves.pdf"
# fig.savefig(img_path)
