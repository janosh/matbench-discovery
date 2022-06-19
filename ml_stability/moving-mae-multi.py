# %%
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from scipy.stats import sem as std_err_of_mean

from ml_stability import ROOT


__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

today = f"{datetime.now():%Y-%m-%d}"

plt.rc("font", size=18)
plt.rc("savefig", bbox="tight", dpi=200)
plt.rcParams["figure.constrained_layout.use"] = True
plt.rc("figure", dpi=150, titlesize=20)


# %%
markers = ["o", "v", "^", "H", "D", ""]

df = pd.read_csv(f"{ROOT}/data/wren-mp-initial-structures.csv").set_index("material_id")


# %%
rare = "all"
# rare = "nla"
# df = df[
#     ~df["composition"].apply(
#         lambda x: any(el.is_rare_earth_metal for el in Composition(x).elements)
#     )
# ]

df_hull = pd.read_csv(f"{ROOT}/data/wbm_e_above_mp_hull.csv").set_index("material_id")

df["e_above_hull"] = df_hull.e_above_hull

df = df.dropna(subset=["e_above_hull"])

tar = df.e_above_hull.to_numpy().ravel()

tar_f = df.filter(like="target").to_numpy().ravel()

pred = df.filter(like="pred").to_numpy().T
mean = np.average(pred, axis=0) - tar_f + tar

res = mean - tar

sort = np.argsort(tar)

tar = tar[sort]
res = res[sort]

half_window = 0.02
increment = 0.002
bot, top = -0.2, 0.3
bins = np.arange(bot, top, increment)

means = np.zeros_like(bins)
std = np.zeros_like(bins)

for j, b in enumerate(bins):
    low = b - half_window
    high = b + half_window

    means[j] = np.mean(np.abs(res[np.argwhere((tar <= high) & (tar > low))]))
    std[j] = std_err_of_mean(np.abs(res[np.argwhere((tar <= high) & (tar > low))]))
fig, ax = plt.subplots(1, figsize=(10, 9))

ax.plot(bins, means)

ax.fill_between(
    bins,
    means + std,
    means - std,
    alpha=0.3,
)

scalebar = AnchoredSizeBar(
    ax.transData,
    2 * half_window,
    "40 meV",
    "lower left",
    pad=0,
    borderpad=0.3,
    # color="white",
    frameon=False,
    size_vertical=0.003,
    # fontproperties=fontprops,
)

ax.add_artist(scalebar)

ax.plot((0.05, 0.5), (0.05, 0.5), color="grey", linestyle="--", alpha=0.3)
ax.plot((-0.5, -0.05), (0.5, 0.05), color="grey", linestyle="--", alpha=0.3)
ax.plot((-0.05, 0.05), (0.05, 0.05), color="grey", linestyle="--", alpha=0.3)
ax.plot((-0.1, 0.1), (0.1, 0.1), color="grey", linestyle="--", alpha=0.3)

ax.fill_between(
    (-0.5, -0.05, 0.05, 0.5),
    (0.5, 0.5, 0.5, 0.5),
    (0.5, 0.05, 0.05, 0.5),
    color="tab:red",
    alpha=0.2,
)

ax.plot((0, 0.05), (0, 0.05), color="grey", linestyle="--", alpha=0.3)
ax.plot((-0.05, 0), (0.05, 0), color="grey", linestyle="--", alpha=0.3)

ax.fill_between(
    (-0.05, 0, 0.05), (0.05, 0.05, 0.05), (0.05, 0, 0.05), color="tab:orange", alpha=0.2
)

ax.annotate(
    xy=(0.055, 0.05),
    xytext=(0.12, 0.05),
    arrowprops=dict(facecolor="black", shrink=0.05),
    text="Corrected\nGGA DFT\nAccuracy",
    verticalalignment="center",
    horizontalalignment="left",
)
ax.annotate(
    xy=(0.105, 0.1),
    xytext=(0.16, 0.1),
    arrowprops=dict(facecolor="black", shrink=0.05),
    text="GGA DFT\nAccuracy",
    verticalalignment="center",
    horizontalalignment="left",
)

ineq = r"|$\Delta E_{Hull-MP}$| > MAE"

ax.text(0, 0.13, ineq, horizontalalignment="center")

ax.set_ylabel("MAE / eV per atom")
ax.set_xlabel(r"$\Delta E_{Hull-MP}$ / eV per atom")

ax.set_ylim((0.0, 0.14))
ax.set_xlim((bot, top))

ax.set_aspect(1.0 / ax.get_data_ratio())


# plt.savefig(f"{PKG_DIR}/plots/{today}-moving-error-wbm-{rare}.pdf")

plt.show()
