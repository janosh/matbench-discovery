# %%
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ml_stability import PKG_DIR, ROOT


__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

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

df["e_above_hull"] = df_hull.e_above_hull


# %%
df = df.dropna(subset=["e_above_hull"])

rare = "all"
target_col = "e_form_target"


pred_mean = df.filter(like="pred").mean(axis=1)
std_epistemic = df.filter(like="pred").var(axis=1, ddof=0)

mean = pred_mean - df[target_col] + df.e_above_hull


std_aleatoric = (df.filter(like="ale") ** 2).mean(axis=1)

std_total = np.sqrt(std_epistemic + std_aleatoric)

# crit = "std"
# test = mean + both

# crit = "neg"
# test = mean - both

crit = "ene"
test = mean

bins = 200
# xlim = (-0.2, 0.2)
xlim = (-0.4, 0.4)
# xlim = (-1, 1)

# thresh = 0.02
thresh = 0.00
# thresh = 0.10
xticks = (-0.4, -0.2, 0, 0.2, 0.4)
# yticks = (0, 300, 600, 900, 1200)

tp = len(df.e_above_hull[(df.e_above_hull <= thresh) & (test <= thresh)])
fn = len(df.e_above_hull[(df.e_above_hull <= thresh) & (test > thresh)])

pos = tp + fn
null = pos / len(df.e_above_hull)

e_type = "true"
tp = df.e_above_hull[(df.e_above_hull <= thresh) & (test <= thresh)]
fn = df.e_above_hull[(df.e_above_hull <= thresh) & (test > thresh)]
fp = df.e_above_hull[(df.e_above_hull > thresh) & (test <= thresh)]
tn = df.e_above_hull[(df.e_above_hull > thresh) & (test > thresh)]
xlabel = r"$\Delta E_{Hull-MP}$ / eV per atom"


# e_type = "pred"
# tp = mean[(tar <= thresh) & (test <= thresh)]
# fn = mean[(tar <= thresh) & (test > thresh)]
# fp = mean[(tar > thresh) & (test <= thresh)]
# tn = mean[(tar > thresh) & (test > thresh)]
# xlabel = r"$\Delta E_{Hull-Pred}$ / eV per atom"
fig, ax = plt.subplots(1, 1, figsize=(10, 9))

ax.hist(
    [tp, fn, fp, tn],
    bins=bins,
    range=xlim,
    alpha=0.5,
    color=["tab:green", "tab:orange", "tab:red", "tab:blue"],
    linestyle="none",
    label=[
        "True Positives",
        "False Negatives",
        "False Positives",
        "True Negatives",
    ],
    stacked=True,
)

ax.legend(frameon=False, loc="upper left")

tp, fp, tn, fn = len(tp), len(fp), len(tn), len(fn)
# null = (tp + fn) / (tp + tn + fp + fn)
ppv = tp / (tp + fp)
tpr = tp / pos
f1 = 2 * ppv * tpr / (ppv + tpr)

print(sum([tp, fp, tn, fn]))

print(f"PPV: {ppv:.2f}")
print(f"TPR: {tpr:.2f}")
print(f"F1: {f1:.2f}")
print(f"Enrich: {ppv/null:.2f}")
print(f"Null: {null:.2f}")

RMSE = ((mean - df[target_col]) ** 2.0).mean() ** 0.5
MAE = (mean - df[target_col]).abs().mean()
print(f"{MAE=:.4}")
print(f"{RMSE=:.4}")

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
# else:
# ax.get_yaxis().set_ticklabels([])

ax.set(xlim=xlim, ylim=ylim)


# NOTE this figure plots hist bars separately which causes aliasing in pdf
# to resolve this take into inkscape and merge regions by colour
plt.savefig(f"{PKG_DIR}/plots/{today}-hist-{e_type}-{crit}-{rare}.pdf")
# plt.savefig(f"{PKG_DIR}/plots/{today}-hist-{e_type}-{crit}-{rare}.png")

plt.show()
