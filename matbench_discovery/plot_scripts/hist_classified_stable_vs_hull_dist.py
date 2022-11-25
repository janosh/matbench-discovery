# %%
import pandas as pd

from matbench_discovery import ROOT, today
from matbench_discovery.plot_scripts import df_wbm
from matbench_discovery.plots import (
    StabilityCriterion,
    WhichEnergy,
    hist_classified_stable_vs_hull_dist,
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


# %%
df = pd.read_csv(
    # f"{ROOT}/data/2022-06-11-from-rhys/wren-mp-initial-structures.csv"
    f"{ROOT}/models/wrenformer/2022-11-15-wrenformer-IS2RE-preds.csv"
).set_index("material_id")

df["e_above_hull"] = df_wbm.e_above_hull_mp2020_corrected_ppd_mp


# %%
nan_counts = df.isna().sum()
assert all(nan_counts == 0), f"df should not have missing values: {nan_counts}"

# target_col = "e_form_target"
target_col = "e_form_per_atom"
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

ax, metrics = hist_classified_stable_vs_hull_dist(
    e_above_hull_pred=df[pred_cols].mean(axis=1) - df[target_col],
    e_above_hull_true=df.e_above_hull,
    which_energy=which_energy,
    stability_crit=stability_crit,
    std_pred=std_total,
    # stability_threshold=-0.05,
    # rolling_accuracy=0,
)

fig = ax.figure
fig.set_size_inches(10, 9)

ax.legend(
    loc="center left",
    frameon=False,
    title=f"Enrichment Factor = {metrics['enrichment']:.3}",
)


# %%
fig_name = f"{today}-wren-wbm-hull-dist-hist-{which_energy=}-{stability_crit=}"
# fig.savefig(f"{ROOT}/figures/{fig_name}.pdf")
