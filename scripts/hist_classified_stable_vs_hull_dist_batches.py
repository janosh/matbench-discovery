# %%
from pymatviz.utils import save_fig

from matbench_discovery import FIGS, today
from matbench_discovery.metrics import (
    df_wbm,
    e_form_col,
    each_pred_col,
    each_true_col,
    stable_metrics,
)
from matbench_discovery.plots import (
    Backend,
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
model_name = "Wrenformer"
which_energy: WhichEnergy = "true"
backend: Backend = "matplotlib"
fig, axs = plt.subplots(2, 3, figsize=(18, 9))
df_wbm[each_pred_col] = df_wbm[each_true_col] + df_wbm[model_name] - df_wbm[e_form_col]


for batch_idx, ax in zip(range(1, 6), axs.flat):
    df_batch = df_wbm[df_wbm.index.str.startswith(f"wbm-{batch_idx}-")]
    assert 1e4 < len(df_batch) < 1e5, print(f"{len(df_batch) = :,}")

    ax = hist_classified_stable_vs_hull_dist(
        df=df_batch,
        each_true_col=each_true_col,
        each_pred_col=each_pred_col,
        which_energy=which_energy,
        ax=ax,
        backend=backend,
    )

    metrics = stable_metrics(df_batch[each_true_col], df_batch[each_pred_col])
    text = f"DAF = {metrics['DAF']:.3}"
    ax.text(0.02, 0.25, text, fontsize=16, transform=ax.transAxes)

    title = f"Batch {batch_idx} ({len(df_batch.filter(like='e_').dropna()):,})"
    ax.set(title=title)


ax = hist_classified_stable_vs_hull_dist(
    df=df_wbm,
    each_true_col=each_true_col,
    each_pred_col=each_pred_col,
    which_energy=which_energy,
    ax=axs.flat[-1],
    backend=backend,
)

metrics = stable_metrics(df_wbm[each_true_col], df_wbm[each_pred_col])
text = f"DAF = {metrics['DAF']:.3}"
ax.text(0.02, 0.3, text, fontsize=16, transform=ax.transAxes)

axs.flat[-1].set(title=f"All batches ({len(df_wbm[model_name].dropna()):,})")
axs.flat[0].legend(frameon=False, loc="upper left")

fig.suptitle(f"{today} {model_name}", y=1.07, fontsize=16)


# %%
img_path = f"{FIGS}/{today}-{model_name}-wbm-hull-dist-hist-batches.pdf"
save_fig(ax, img_path)
