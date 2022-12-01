# %%
from matbench_discovery import ROOT, today
from matbench_discovery.load_preds import load_df_wbm_with_preds
from matbench_discovery.plots import (
    WhichEnergy,
    hist_classified_stable_vs_hull_dist,
    plt,
)

__author__ = "Janosh Riebesell"
__date__ = "2022-12-01"

"""
Histogram of the energy difference (either according to DFT ground truth [default] or
model predicted energy) to the convex hull for materials in the WBM data set. The
histogram stacks true/false positives/negatives with different colors.
"""


# %%
models = (
    "Wren, CGCNN, CGCNN IS2RE, CGCNN RS2RE, Voronoi RF, "
    "Wrenformer, MEGNet, M3GNet, BOWSR MEGNet"
).split(", ")
df_wbm = load_df_wbm_with_preds(models=models).round(3)

target_col = "e_form_per_atom_mp2020_corrected"
e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"


# %%
which_energy: WhichEnergy = "true"
fig, axs = plt.subplots(3, 3, figsize=(18, 12))

model_name = "Wrenformer"

for model_name, ax in zip(models, axs.flat, strict=True):

    ax, metrics = hist_classified_stable_vs_hull_dist(
        e_above_hull_true=df_wbm[e_above_hull_col],
        e_above_hull_pred=df_wbm[e_above_hull_col]
        + (df_wbm[model_name] - df_wbm[target_col]),
        which_energy=which_energy,
        ax=ax,
    )

    text = f"Enrichment\nFactor = {metrics['enrichment']:.3}"
    ax.text(0.02, 0.25, text, fontsize=16, transform=ax.transAxes)

    title = f"{model_name} ({len(df_wbm[model_name].dropna()):,})"
    ax.set(title=title)


# axs.flat[0].legend(frameon=False, loc="upper left")

fig.suptitle(f"{today} {which_energy=}", y=1.07, fontsize=16)


# %%
img_path = f"{ROOT}/figures/{today}-wbm-hull-dist-hist-models.pdf"
ax.figure.savefig(img_path)
