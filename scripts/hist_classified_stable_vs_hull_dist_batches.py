"""Histogram of the energy difference (either according to DFT ground truth [default] or
model predicted energy) to the convex hull for materials in the WBM data set. The
histogram stacks true/false positives/negatives with different colors.

See fig. S1 in https://science.org/doi/10.1126/sciadv.abn4117.
"""

# %%
from typing import Final

import pandas as pd
import pymatviz as pmv
from pymatviz.enums import Key
from pymatviz.utils import MATPLOTLIB, PLOTLY

from matbench_discovery import PDF_FIGS
from matbench_discovery.data import Model
from matbench_discovery.enums import MbdKey
from matbench_discovery.plots import hist_classified_stable_vs_hull_dist
from matbench_discovery.preds import df_preds

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-08-25"


# %%
model_name = Model.mace.label
which_energy: Final = "true"
backend: Final = MATPLOTLIB
df_preds[Key.each_pred] = (
    df_preds[MbdKey.each_true] + df_preds[model_name] - df_preds[MbdKey.e_form_dft]
)
df_preds[(batch_col := "batch_idx")] = df_preds.index.str.split("-").str[1].astype(int)


# %% matplotlib
# fig, axs = plt.subplots(2, 3, figsize=(18, 9))
# for batch_idx, ax in zip(range(1, 6), axs.flat):
#     df_batch = df_wbm[df_wbm.index.str.startswith(f"wbm-{batch_idx}-")]
#     assert 1e4 < len(df_batch) < 1e5, print(f"{len(df_batch) = :,}")

#     ax = hist_classified_stable_vs_hull_dist(
#         df=df_batch,
#         each_true_col=MbdKey.each_true,
#         each_pred_col=Key.each_pred,
#         which_energy=which_energy,
#         ax=ax,
#         backend=backend,
#     )

#     metrics = stable_metrics(df_batch[MbdKey.each_true], df_batch[Key.each_pred])
#     text = f"DAF = {metrics['DAF']:.3}"
#     ax.text(0.02, 0.25, text, fontsize=16, transform=ax.transAxes)

#     title = f"Batch {batch_idx} ({len(df_batch.filter(like='e_').dropna()):,})"
#     ax.set(title=title)


# ax = hist_classified_stable_vs_hull_dist(
#     df=df_wbm,
#     each_true_col=MbdKey.each_true,
#     each_pred_col=Key.each_pred,
#     which_energy=which_energy,
#     ax=axs.flat[-1],
#     backend=backend,
# )

# metrics = stable_metrics(df_wbm[MbdKey.each_true], df_wbm[Key.each_pred])
# text = f"DAF = {metrics['DAF']:.3}"
# ax.text(0.02, 0.3, text, fontsize=16, transform=ax.transAxes)

# axs.flat[-1].set(title=f"All batches ({len(df_wbm[model_name].dropna()):,})")
# axs.flat[0].legend(frameon=False, loc="upper left")

# fig.suptitle(f"{today} {model_name}", y=1.07, fontsize=16)


# %% plotly
fig = hist_classified_stable_vs_hull_dist(
    # # plot whole df as last subplot
    df=pd.concat([df_preds, df_preds.assign(batch_idx=6)]),
    each_true_col=MbdKey.each_true,
    each_pred_col=Key.each_pred,
    which_energy=which_energy,
    facet_col=batch_col,
    facet_col_wrap=2,
    stability_threshold=None,
    rolling_acc=None,
    backend=PLOTLY,
)
for anno in fig.layout.annotations:
    if not anno.text.startswith("batch_idx="):
        continue
    batch_idx = int(anno.text.split("=", 1)[-1])
    len_df = sum(df_preds[batch_col] == int(batch_idx))
    anno.text = f"Batch {batch_idx} ({len_df:,})"

# allow scrolling and zooming each subplot individually
fig.update_xaxes(matches=None, range=[-0.5, 0.5])
# remove y-axis labels
fig.update_yaxes(matches=None, title="")
fig.show()


# %%
img_path = f"{PDF_FIGS}/{model_name}-wbm-hull-dist-hist-batches.pdf"
pmv.save_fig(fig, img_path)
