# %%
from datetime import datetime

import matplotlib.pyplot as plt
import pandas as pd

from ml_stability import PKG_DIR, ROOT


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

plt.rc("font", size=18)
plt.rc("savefig", bbox="tight", dpi=200)
plt.rcParams["figure.constrained_layout.use"] = True
plt.rc("figure", dpi=150, titlesize=20)


# %%
df = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wren-mp-initial-structures.csv"
).set_index("material_id")

df_hull = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")

df["e_above_mp_hull"] = df_hull.e_above_mp_hull


df_summary = pd.read_csv(f"{ROOT}/data/wbm-steps-summary.csv", comment="#").set_index(
    "material_id"
)


# %%
assert df.e_above_mp_hull.isna().sum() == 0

rare = "all"
target_col = "e_form_target"

var_aleatoric = (df.filter(like="ale") ** 2).mean(axis=1)
var_epistemic = df.filter(like="pred").var(axis=1, ddof=0)

std_total = (var_epistemic + var_aleatoric) ** 0.5

criterion = "energy"
error = df.filter(like="pred").mean(axis=1) - df[target_col]
mean = error + df.e_above_mp_hull

test = mean
# criterion = "std"
# test += std_total

# criterion = "neg"
# test -= std_total

xlim = (-0.4, 0.4)

# set stability threshold at on or 0.1 eV / atom above the hull
stability_thresh = (0, 0.1)[0]
xticks = (-0.4, -0.2, 0, 0.2, 0.4)
# yticks = (0, 300, 600, 900, 1200)

# --- histogram by DFT-computed distance to convex hull
e_type = "true"
actual_pos = df.e_above_mp_hull <= stability_thresh
actual_neg = df.e_above_mp_hull > stability_thresh
model_pos = test <= stability_thresh
model_neg = test > stability_thresh

n_true_pos = len(df.e_above_mp_hull[actual_pos & model_pos])
n_false_neg = len(df.e_above_mp_hull[actual_pos & model_neg])

n_total_pos = n_true_pos + n_false_neg
null = n_total_pos / len(df.e_above_mp_hull)

true_pos = df.e_above_mp_hull[actual_pos & model_pos]
false_neg = df.e_above_mp_hull[actual_pos & model_neg]
false_pos = df.e_above_mp_hull[actual_neg & model_pos]
true_neg = df.e_above_mp_hull[actual_neg & model_neg]
xlabel = r"$\Delta E_{Hull-MP}$ / eV per atom"


# --- histogram by model-predicted distance to convex hull
# e_type = "pred"
# true_pos = mean[actual_pos & model_pos]
# false_neg = mean[actual_pos & model_neg]
# false_pos = mean[actual_neg & model_pos]
# true_neg = mean[actual_neg & model_neg]
# xlabel = r"$\Delta E_{Hull-Pred}$ / eV per atom"
fig, ax = plt.subplots(1, 1, figsize=(10, 9))

ax.hist(
    [true_pos, false_neg, false_pos, true_neg],
    bins=200,
    range=xlim,
    alpha=0.5,
    color=["tab:green", "tab:orange", "tab:red", "tab:blue"],
    label=["True Positives", "False Negatives", "False Positives", "True Negatives"],
    stacked=True,
)

ax.legend(frameon=False, loc="upper left")

n_true_pos, n_false_pos, n_true_neg, n_false_neg = (
    len(true_pos),
    len(false_pos),
    len(true_neg),
    len(false_neg),
)
# null = (tp + fn) / (tp + tn + fp + fn)
ppv = n_true_pos / (n_true_pos + n_false_pos)
tpr = n_true_pos / n_total_pos
f1 = 2 * ppv * tpr / (ppv + tpr)

assert n_true_pos + n_false_pos + n_true_neg + n_false_neg == len(df)

print(f"PPV: {ppv:.2f}")
print(f"TPR: {tpr:.2f}")
print(f"F1: {f1:.2f}")
print(f"Enrich: {ppv/null:.2f}")
print(f"Null: {null:.2f}")

RMSE = (error**2.0).mean() ** 0.5
MAE = error.abs().mean()
print(f"{MAE=:.3}")
print(f"{RMSE=:.3}")

ylim = (0, 6000)
yticks = (0, 2000, 4000, 6000)

xpos, ypos = 0.45 * xlim[1], 0.96 * ylim[1]
fontsize = 20

ax.text(
    xpos,
    ypos,
    # f"Prevalence = {null:.2f}\nPrecision = {ppv:.2f}\nRecall = {tpr:.2f}",
    f"Enrichment\nFactor = {ppv/null:.1f}",
    fontsize=fontsize,
    verticalalignment="top",
    horizontalalignment="left",
)

xpos, ypos = 0.90 * xlim[0], 0.96 * ylim[1]


ax.set(xticks=xticks, yticks=yticks)
ax.set(xlabel=xlabel, ylabel="Number of Compounds")
# ax.get_yaxis().set_ticklabels([])

ax.set(xlim=xlim, ylim=ylim)


# NOTE this figure plots hist bars separately which causes aliasing in pdf
# to resolve this take into Inkscape and merge regions by color
img_path = f"{PKG_DIR}/plots/{today}-hist-{e_type=}-{criterion=}-{rare=}.pdf"
# plt.savefig(img_path)
