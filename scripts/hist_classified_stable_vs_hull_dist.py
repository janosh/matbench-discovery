# %%
from pymatviz.utils import save_fig

from matbench_discovery import FIGS, today
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.plots import WhichEnergy, hist_classified_stable_vs_hull_dist

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

"""
Histogram of the energy difference (either according to DFT ground truth [default] or
model predicted energy) to the convex hull for materials in the WBM data set. The
histogram stacks true/false positives/negatives with different colors.

See fig. S1 in https://science.org/doi/10.1126/sciadv.abn4117.
"""


# %%
model_name = "Wrenformer"
df_wbm = load_df_wbm_with_preds(models=[model_name]).round(3)


# %%
target_col = "e_form_per_atom_mp2020_corrected"
e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"
which_energy: WhichEnergy = "true"
# std_factor=0,+/-1,+/-2,... changes the criterion for material stability to
# energy+std_factor*std. energy+std means predicted energy plus the model's uncertainty
# in the prediction have to be on or below the convex hull to be considered stable. This
# reduces the false positive rate, but increases the false negative rate. Vice versa for
# energy-std. energy+std should be used for cautious exploration, energy-std for
# exhaustive exploration.
std_factor = 0

# TODO column names to compute standard deviation from are currently hardcoded
# needs to be updated when adding non-aviary models with uncertainty estimation
var_aleatoric = (df_wbm.filter(like="_ale_") ** 2).mean(axis=1)
var_epistemic = df_wbm.filter(regex=r"_pred_\d").var(axis=1, ddof=0)
std_total = (var_epistemic + var_aleatoric) ** 0.5
std_total = df_wbm[f"{model_name}_std"]
e_above_hull_pred = (
    df_wbm[e_above_hull_col]
    + (df_wbm[model_name] + std_factor * std_total)
    - df_wbm[target_col]
)

ax, metrics = hist_classified_stable_vs_hull_dist(
    e_above_hull_true=df_wbm[e_above_hull_col],
    e_above_hull_pred=e_above_hull_pred,
    which_energy=which_energy,
    # stability_threshold=-0.05,
    rolling_accuracy=0,
    # backend="matplotlib",
)
if hasattr(ax, "legend"):
    legend_title = f"Enrichment Factor = {metrics['enrichment']:.3}"
    ax.legend(loc="upper left", frameon=False, title=legend_title)

ax


# %%
img_path = f"{FIGS}/{today}-wren-wbm-hull-dist-hist-{which_energy=}"
# save_fig(ax, f"{img_path}.pdf")
save_fig(ax, f"{img_path}.html")
