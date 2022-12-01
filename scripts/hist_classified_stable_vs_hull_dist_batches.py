# %%
from matbench_discovery import ROOT, today
from matbench_discovery.load_preds import load_df_wbm_with_preds
from matbench_discovery.plots import (
    WhichEnergy,
    hist_classified_stable_vs_hull_dist,
    plt,
)

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-08-25"

"""
Histogram of the energy difference (either according to DFT ground truth [default] or
model predicted energy) to the convex hull for materials in the WBM data set. The
histogram stacks true/false positives/negatives with different colors.

See fig. S1 in https://science.org/doi/10.1126/sciadv.abn4117.
"""


# %%
df_wbm = load_df_wbm_with_preds(models="Wren Wrenformer".split()).round(3)
target_col = "e_form_per_atom_mp2020_corrected"
e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"


# %%
which_energy: WhichEnergy = "true"
fig, axs = plt.subplots(2, 3, figsize=(18, 9))

model_name = "Wrenformer"

for batch_idx, ax in zip(range(1, 6), axs.flat):
    batch_df = df_wbm[df_wbm.index.str.startswith(f"wbm-step-{batch_idx}-")]
    assert 1e4 < len(batch_df) < 1e5, print(f"{len(batch_df) = :,}")

    ax, metrics = hist_classified_stable_vs_hull_dist(
        e_above_hull_true=batch_df[e_above_hull_col],
        e_above_hull_pred=batch_df[e_above_hull_col]
        + (batch_df[model_name] - batch_df[target_col]),
        which_energy=which_energy,
        ax=ax,
    )

    text = f"Enrichment\nFactor = {metrics['enrichment']:.3}"
    ax.text(0.02, 0.25, text, fontsize=16, transform=ax.transAxes)

    title = f"Batch {batch_idx} ({len(batch_df.filter(like='e_').dropna()):,})"
    ax.set(title=title)


ax, metrics = hist_classified_stable_vs_hull_dist(
    e_above_hull_true=df_wbm[e_above_hull_col],
    e_above_hull_pred=df_wbm[e_above_hull_col]
    + (df_wbm[model_name] - df_wbm[target_col]),
    which_energy=which_energy,
    ax=axs.flat[-1],
)

text = f"Enrichment\nFactor = {metrics['enrichment']:.3}"
ax.text(0.02, 0.3, text, fontsize=16, transform=ax.transAxes)

axs.flat[-1].set(title=f"All batches ({len(df_wbm[model_name].dropna()):,})")
axs.flat[0].legend(frameon=False, loc="upper left")

fig.suptitle(f"{today} {model_name}", y=1.07, fontsize=16)


# %%
img_path = f"{ROOT}/figures/{today}-{model_name}-wbm-hull-dist-hist-batches.pdf"
# ax.figure.savefig(img_path)
