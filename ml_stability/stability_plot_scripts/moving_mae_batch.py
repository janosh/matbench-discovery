# %%
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

from ml_stability import PKG_DIR, ROOT


__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

today = f"{datetime.now():%Y-%m-%d}"


plt.rc("font", size=18)
plt.rc("savefig", bbox="tight", dpi=200)
plt.rcParams["figure.constrained_layout.use"] = True
plt.rc("figure", dpi=150, titlesize=20)


# %%
rare = "all"

df_wbm = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wren-mp-initial-structures.csv"
).set_index("material_id")

df_hull = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")

df_wbm["e_above_hull"] = df_hull.e_above_hull


# %%
fig, ax = plt.subplots(1, figsize=(10, 9))

markers = [
    "o",
    "v",
    "^",
    "H",
    "D",
    # "",
]
df_wbm = df_wbm.dropna(subset=["e_above_hull"])

for i, m in enumerate(markers):
    offsets = 1
    title = f"Batch-{i+offsets}"

    df = df_wbm[df_wbm.index.str.contains(f"wbm-step-{i+offsets}")]
    tar = df.e_above_hull.to_numpy().ravel()

    tar_f = df.filter(like="target").to_numpy().ravel()

    pred = df.filter(like="pred").to_numpy().T
    # mean = np.average(pred, axis=0)
    mean = np.average(pred, axis=0) - tar_f + tar

    res = np.abs(mean - tar)

    # sort = np.argsort(tar)

    # tar = tar[sort]
    # res = res[sort]

    # tar = mean

    # half_window = 0.1
    # half_window = 0.01
    half_window = 0.02
    increment = 0.002
    bottom, top = -0.2, 0.3
    # bot, top = -0.2, 0.6
    bins = np.arange(bottom, top, increment)

    means = np.zeros_like(bins)
    medians = np.zeros_like(bins)
    quant = np.zeros_like(bins)

    for j, b in enumerate(bins):
        low = b - half_window
        high = b + half_window

        means[j] = np.mean(res[np.argwhere((tar <= high) & (tar > low))])
        # medians[j] = np.median(res[np.argwhere((tar <= high) & (tar > low))])
        # quant[j] = np.nanquantile(res[np.argwhere((tar <= high) & (tar > low))], 0.9)

    ax.plot(
        bins,
        means,
        label=title,
        marker=m,
        markevery=20,
        markersize=8 if m == "D" else 10,
        # fillstyle="none",
        fillstyle="full",
        markerfacecolor="white",
        mew=2.5,
    )
    # ax.plot(bins, medians, label=title)
    # ax.plot(bins, quant, label=title)


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

# ax.fill_between(
#     (-0.5, -0.05, 0.05, 0.5),
#     (0.5, 0.5, 0.5, 0.5),
#     (0.5, 0.05, 0.05, 0.5),
#     color="tab:red",
#     alpha=0.2,
# )

ax.plot((0, 0.05), (0, 0.05), color="grey", linestyle="--", alpha=0.3)
ax.plot((-0.05, 0), (0.05, 0), color="grey", linestyle="--", alpha=0.3)

# ax.fill_between(
#     (-0.05, 0, 0.05),
#     (0.05, 0.05, 0.05),
#     (0.05, 0, 0.05),
#     color="tab:orange",
#     alpha=0.2,
# )

# ax.fill_between((-0.5, 0), (0.5, 0), (0, 0), color="tab:green", alpha=0.2)
# ax.fill_between((0, 0.5), (0, 0.5), (0, 0), color="tab:green", alpha=0.2)

ax.set(xlabel=r"$\Delta E_{Hull-MP}$ / eV per atom", ylabel="MAE / eV per atom")

ax.set(xlim=(bottom, top), ylim=(0.0, 0.14))
ax.legend(
    # title=r"$\bf{Wren}$",
    # frameon=False,
    loc="lower right",
    # loc="upper left",
    facecolor="white",
    framealpha=1.0,
    edgecolor="white",
)

# titles = [r"Wren\ (This\ Work)", "CGCNN-D"]
# reps = [
#     r"Wyckoff\ Representation",
#     "Pre \u2212 relaxation\\ Structures",
# ]

# ax.annotate(
#     f"$\\bf{{Input: {reps[0]}}}$\n$\\bf{{Model: {titles[0]}}}$\n ",
#     (0.05, 0.82),
#     xycoords="axes fraction",
# )


plt.savefig(f"{PKG_DIR}/plots/{today}-moving-error-wbm-{rare}-batches.pdf")
# plt.savefig(f"{PKG_DIR}/plots/{today}-moving-error-wbm-{rare}-batches.pdf")


plt.show()
