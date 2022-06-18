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
fig, ax = plt.subplots(1, 1, figsize=(10, 9))


df = pd.read_csv(f"{ROOT}/data/wren-mp-initial-structures.csv")

df_hull = pd.read_csv(f"{ROOT}/data/wbm_e_above_mp.csv")

df["e_above_hull"] = pd.to_numeric(
    df["material_id"].map(dict(zip(df_hull.material_id, df_hull.e_above_hull)))
)


# %%
df = df.dropna(axis=0, subset=["e_above_hull"])

init = len(df)


rare = "all"

# rare = "nla"
# df = df[
#     ~df["composition"].apply(
#         lambda x: any(el.is_rare_earth_metal for el in Composition(x).elements)
#     )
# ]

# print(1-len(df)/init)

tar = df.e_above_hull.to_numpy().ravel()

print(len(tar))

# tar = df.filter(like="target").to_numpy().ravel() - e_hull
tar_f = df.filter(like="target").to_numpy().ravel()

pred = df.filter(like="pred").to_numpy().T
# mean = np.average(pred, axis=0) - e_hull
mean = np.average(pred, axis=0) - tar_f + tar

epi = np.var(pred, axis=0, ddof=0)

ales = df.filter(like="ale").to_numpy().T
ale = np.mean(np.square(ales), axis=0)

both = np.sqrt(epi + ale)

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

tp = len(tar[(tar <= thresh) & (test <= thresh)])
fn = len(tar[(tar <= thresh) & (test > thresh)])

pos = tp + fn
null = pos / len(tar)

e_type = "true"
tp = tar[(tar <= thresh) & (test <= thresh)]
fn = tar[(tar <= thresh) & (test > thresh)]
fp = tar[(tar > thresh) & (test <= thresh)]
tn = tar[(tar > thresh) & (test > thresh)]
xlabel = r"$\Delta E_{Hull-MP}$ / eV per atom"


# e_type = "pred"
# tp = mean[(tar <= thresh) & (test <= thresh)]
# fn = mean[(tar <= thresh) & (test > thresh)]
# fp = mean[(tar > thresh) & (test <= thresh)]
# tn = mean[(tar > thresh) & (test > thresh)]
# xlabel = r"$\Delta E_{Hull-Pred}$ / eV per atom"

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

print(f"PPV: {ppv}")
print(f"TPR: {tpr}")
print(f"F1: {f1}")
print(f"Enrich: {ppv/null}")
print(f"Null: {null}")

# print(f"MAE: {np.mean(np.abs(mean - tar))}")
# print(f"RMSE: {np.sqrt(np.mean(np.square(mean - tar)))}")

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


ax.set_xticks(xticks)
ax.set_yticks(yticks)

ax.set_xlabel(xlabel)

ax.set_ylabel("Number of Compounds")
# else:
# ax.get_yaxis().set_ticklabels([])

ax.set_ylim(ylim)
ax.set_xlim(xlim)

ax.set_aspect(1.0 / ax.get_data_ratio())


# NOTE this figure plots hist bars separately which causes aliasing in pdf
# to resolve this take into inkscape and merge regions by colour
plt.savefig(f"{PKG_DIR}/plots/{today}-hist-{e_type}-{crit}-{rare}.pdf")
# plt.savefig(f"{PKG_DIR}/plots/{today}-hist-{e_type}-{crit}-{rare}.png")

plt.show()
