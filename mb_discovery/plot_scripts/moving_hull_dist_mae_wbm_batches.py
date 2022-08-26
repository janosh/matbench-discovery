# %%
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

from mb_discovery import ROOT


__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

today = f"{datetime.now():%Y-%m-%d}"


plt.rc("savefig", bbox="tight", dpi=200)
plt.rcParams["figure.constrained_layout.use"] = True
plt.rc("figure", dpi=150)
plt.rc("font", size=16)


# %%
rare = "all"

df_wbm = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wren-mp-initial-structures.csv"
).set_index("material_id")

df_hull = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")

df_wbm["e_above_mp_hull"] = df_hull.e_above_mp_hull
assert df_wbm.e_above_mp_hull.isna().sum() == 0

target_col = "e_form_target"

# make sure we average the expected number of ensemble member predictions
assert df_wbm.filter(regex=r"_pred_\d").shape[1] == 10

df_wbm["e_above_mp_hull_pred"] = (
    df_wbm.filter(regex=r"_pred_\d").mean(axis=1)
    - df_wbm[target_col]
    + df_wbm.e_above_mp_hull
)
df_wbm["error"] = abs(df_wbm.e_above_mp_hull_pred - df_wbm.e_above_mp_hull)


# %%
fig, ax = plt.subplots(1, figsize=(10, 9))
markers = ("o", "v", "^", "H", "D")
assert len(markers) == 5  # number of WBM round of element substitution

for idx, marker in enumerate(markers):
    title = f"Batch-{idx}"

    df = df_wbm[df_wbm.index.str.startswith(f"wbm-step-{idx}")]

    # half_window = 0.1
    # half_window = 0.01
    half_window = 0.02
    increment = 0.002
    bottom, top = -0.2, 0.3
    # bottom, top = -0.2, 0.6
    bins = np.arange(bottom, top, increment)

    means = np.zeros_like(bins)
    medians = np.zeros_like(bins)
    quant = np.zeros_like(bins)

    for jdx, bin_center in enumerate(bins):
        low = bin_center - half_window
        high = bin_center + half_window

        mask = (df.e_above_mp_hull <= high) & (df.e_above_mp_hull > low)

        means[jdx] = df.error[mask].mean()
        # medians[jdx] = df.error[mask].median()
        # quant[jdx] = np.nanquantile(df.error[mask], 0.9)

    ax.plot(
        bins,
        means,
        label=title,
        marker=marker,
        markevery=20,
        markersize=8 if marker == "D" else 10,
        # fillstyle="none",
        fillstyle="full",
        markerfacecolor="white",
        mew=2.5,
    )
    # ax.plot(bins, medians, label=title)
    # ax.plot(bins, quant, label=title)


scale_bar = AnchoredSizeBar(
    ax.transData,
    2 * half_window,
    "40 meV",
    "lower left",
    pad=1,
    borderpad=0.3,
    # color="white",
    frameon=False,
    size_vertical=0.003,
)

ax.add_artist(scale_bar)

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

img_path = f"{ROOT}/figures/{today}-moving-hull-dist-mae-wbm-batches-{rare=}.pdf"
# plt.savefig(img_path)


plt.show()
