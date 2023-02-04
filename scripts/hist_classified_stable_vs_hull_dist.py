# %%
from pymatviz.utils import save_fig

from matbench_discovery import FIGS, today
from matbench_discovery.metrics import stable_metrics
from matbench_discovery.plots import WhichEnergy, hist_classified_stable_vs_hull_dist
from matbench_discovery.preds import df_wbm, e_form_col, each_pred_col, each_true_col

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
df_wbm[each_pred_col] = df_wbm[each_true_col] + (
    (df_wbm[model_name] + std_factor * std_total) - df_wbm[e_form_col]
)

fig = hist_classified_stable_vs_hull_dist(
    df_wbm,
    each_true_col=each_true_col,
    each_pred_col=each_pred_col,
    which_energy=which_energy,
    # stability_threshold=-0.05,
    # rolling_acc=0,
    backend="plotly",
)

metrics = stable_metrics(df_wbm[each_true_col], df_wbm[each_pred_col])
legend_title = f"DAF = {metrics['DAF']:.3}"

if hasattr(fig, "legend"):  # matplotlib
    fig.legend(loc="upper left", frameon=False, title=legend_title)
else:  # plotly
    fig.layout.legend.title.text = legend_title
    fig.show()


# %%
img_path = f"{FIGS}/{today}-wren-wbm-hull-dist-hist-{which_energy=}"
# save_fig(fig, f"{img_path}.svelte")
save_fig(fig, f"{img_path}.webp")
