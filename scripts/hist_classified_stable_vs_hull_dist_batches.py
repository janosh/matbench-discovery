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

from matbench_discovery import PDF_FIGS
from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.plots import hist_classified_stable_vs_hull_dist
from matbench_discovery.preds.discovery import df_preds

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-08-25"


# %%
model_name = Model.mace_mp_0.label
which_energy: Final = "true"
df_preds[Key.each_pred] = (
    df_preds[MbdKey.each_true] + df_preds[model_name] - df_preds[MbdKey.e_form_dft]
)
df_preds[(batch_col := "batch_idx")] = df_preds.index.str.split("-").str[1].astype(int)


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
