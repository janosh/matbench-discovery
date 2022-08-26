# %%
from datetime import datetime

import matplotlib.pyplot as plt
import pandas as pd

from mb_discovery import ROOT
from mb_discovery.plot_scripts.plot_funcs import (
    hist_classified_stable_as_func_of_hull_dist,
)


__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-08-25"

"""
Histogram of the energy difference (either according to DFT ground truth [default] or
model predicted energy) to the convex hull for materials in the WBM data set. The
histogram is broken down into true positives, false negatives, false positives, and true
negatives based on whether the model predicts candidates to be below the known convex
hull. Ideally, in discovery setting a model should exhibit high recall, i.e. the
majority of materials below the convex hull being correctly identified by the model.

See fig. S1 in https://science.org/doi/10.1126/sciadv.abn4117.
"""

today = f"{datetime.now():%Y-%m-%d}"

plt.rc("savefig", bbox="tight", dpi=200)
plt.rcParams["figure.constrained_layout.use"] = True
plt.rc("figure", dpi=150)
plt.rc("font", size=16)


# %%
df = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wren-mp-initial-structures.csv"
).set_index("material_id")

df_hull = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")

df["e_above_mp_hull"] = df_hull.e_above_mp_hull

# download wbm-steps-summary.csv (23.31 MB)
df_summary = pd.read_csv(
    "https://figshare.com/ndownloader/files/36714216?private_link=ff0ad14505f9624f0c05"
).set_index("material_id")


# %%
assert df.e_above_mp_hull.isna().sum() == 0

energy_type = "true"
criterion = "energy"
df["wbm_batch"] = df.index.str.split("-").str[2]
fig, axs = plt.subplots(2, 3, figsize=(18, 9))

# make sure we average the expected number of ensemble member predictions
pred_cols = df.filter(regex=r"_pred_\d").columns
assert len(pred_cols) == 10

common_kwargs = dict(
    target_col="e_form_target",
    pred_cols=pred_cols,
    energy_type=energy_type,
    criterion=criterion,
    e_above_hull_col="e_above_mp_hull",
)

for (batch_idx, batch_df), ax in zip(df.groupby("wbm_batch"), axs.flat):
    hist_classified_stable_as_func_of_hull_dist(batch_df, ax=ax, **common_kwargs)

    title = f"Batch {batch_idx} ({len(df):,})"
    ax.set(title=title)


hist_classified_stable_as_func_of_hull_dist(df, ax=axs.flat[-1], **common_kwargs)

axs.flat[-1].set(title=f"Combined {batch_idx} ({len(df):,})")
axs.flat[0].legend(frameon=False, loc="upper left")

img_name = f"{today}-wren-wbm-hull-dist-hist-{energy_type=}-{criterion=}.pdf"
# plt.savefig(f"{ROOT}/figures/{img_name}")
