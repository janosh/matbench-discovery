# %%
from datetime import datetime
from typing import Literal

import matplotlib.pyplot as plt
import pandas as pd

from mb_discovery import ROOT
from mb_discovery.plot_scripts.plot_funcs import (
    hist_classified_stable_as_func_of_hull_dist,
)


__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

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
plt.rc("figure", dpi=200)
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
nan_counts = df.isna().sum()
assert all(nan_counts == 0), f"df should not have missing values: {nan_counts}"

target_col = "e_form_target"
stability_crit: Literal["energy", "energy+std", "energy-std"] = "energy"
energy_type: Literal["true", "pred"] = "true"


# make sure we average the expected number of ensemble member predictions
pred_cols = df.filter(regex=r"_pred_\d").columns
assert len(pred_cols) == 10

ax = hist_classified_stable_as_func_of_hull_dist(
    df,
    target_col,
    pred_cols,
    e_above_hull_col="e_above_mp_hull",
    energy_type=energy_type,
    stability_crit=stability_crit,
)

ax.figure.set_size_inches(10, 9)

ax.legend(loc="upper left", frameon=False)

fig_name = f"wren-wbm-hull-dist-hist-{energy_type=}-{stability_crit=}"
img_path = f"{ROOT}/figures/{today}-{fig_name}.pdf"
# plt.savefig(img_path)
