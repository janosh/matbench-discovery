# %%
from datetime import datetime

import matplotlib.pyplot as plt
import pandas as pd

from mb_discovery import ROOT
from mb_discovery.plot_scripts.plot_funcs import (
    StabilityCriterion,
    WhichEnergy,
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
    "https://figshare.com/files/36714216?private_link=ff0ad14505f9624f0c05"
).set_index("material_id")


# %%
nan_counts = df.isna().sum()
assert all(nan_counts == 0), f"df should not have missing values: {nan_counts}"

target_col = "e_form_target"
stability_crit: StabilityCriterion = "energy"
which_energy: WhichEnergy = "true"

if "std" in stability_crit:
    # TODO column names to compute standard deviation from are currently hardcoded
    # needs to be updated when adding non-aviary models with uncertainty estimation
    var_aleatoric = (df.filter(like="_ale_") ** 2).mean(axis=1)
    var_epistemic = df.filter(regex=r"_pred_\d").var(axis=1, ddof=0)
    std_total = (var_epistemic + var_aleatoric) ** 0.5
else:
    std_total = None

# make sure we average the expected number of ensemble member predictions
pred_cols = df.filter(regex=r"_pred_\d").columns
assert len(pred_cols) == 10

ax = hist_classified_stable_as_func_of_hull_dist(
    e_above_hull_pred=df[pred_cols].mean(axis=1) - df[target_col],
    e_above_hull_true=df.e_above_mp_hull,
    which_energy=which_energy,
    stability_crit=stability_crit,
    std_pred=std_total,
)

ax.figure.set_size_inches(10, 9)

ax.legend(loc="upper left", frameon=False)

fig_name = f"wren-wbm-hull-dist-hist-{which_energy=}-{stability_crit=}"
img_path = f"{ROOT}/figures/{today}-{fig_name}.pdf"
# plt.savefig(img_path)
