# %%
from datetime import datetime

import pandas as pd

from mb_discovery import ROOT
from mb_discovery.plots import (
    StabilityCriterion,
    WhichEnergy,
    hist_classified_stable_as_func_of_hull_dist,
    plt,
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
    "https://figshare.com/files/37570234?private_link=ff0ad14505f9624f0c05"
).set_index("material_id")


# %%
which_energy: WhichEnergy = "true"
stability_crit: StabilityCriterion = "energy"
df["wbm_batch"] = df.index.str.split("-").str[2]
fig, axs = plt.subplots(2, 3, figsize=(18, 9))

# make sure we average the expected number of ensemble member predictions
pred_cols = df.filter(regex=r"_pred_\d").columns
assert len(pred_cols) == 10


for (batch_idx, batch_df), ax in zip(df.groupby("wbm_batch"), axs.flat):
    hist_classified_stable_as_func_of_hull_dist(
        e_above_hull_pred=batch_df[pred_cols].mean(axis=1) - batch_df.e_form_target,
        e_above_hull_true=batch_df.e_above_mp_hull,
        which_energy=which_energy,
        stability_crit=stability_crit,
        ax=ax,
    )

    title = f"Batch {batch_idx} ({len(df):,})"
    ax.set(title=title)


hist_classified_stable_as_func_of_hull_dist(
    e_above_hull_pred=df[pred_cols].mean(axis=1),
    e_above_hull_true=df.e_above_mp_hull,
    which_energy=which_energy,
    stability_crit=stability_crit,
    ax=axs.flat[-1],
)

axs.flat[-1].set(title=f"Combined {batch_idx} ({len(df):,})")
axs.flat[0].legend(frameon=False, loc="upper left")

img_name = f"{today}-wren-wbm-hull-dist-hist-{which_energy=}-{stability_crit=}.pdf"
# plt.savefig(f"{ROOT}/figures/{img_name}")
