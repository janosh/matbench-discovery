# %%
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from scipy.stats import sem as std_err_of_mean

from ml_stability import PKG_DIR, ROOT


__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

today = f"{datetime.now():%Y-%m-%d}"

plt.rc("font", size=18)
plt.rc("savefig", bbox="tight", dpi=200)
plt.rcParams["figure.constrained_layout.use"] = True
plt.rc("figure", dpi=150, titlesize=20)


# %%
markers = ["o", "v", "^", "H", "D", ""]

df = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wren-mp-initial-structures.csv"
).set_index("material_id")


# %%
rare = "all"
# from pymatgen.core import Composition
# rare = "no-lanthanides"
# df["contains_rare_earths"] = df.composition.map(
#     lambda x: any(el.is_rare_earth_metal for el in Composition(x))
# )
# df = df.query("~contains_rare_earths")


df_hull = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")

df["e_above_mp_hull"] = df_hull.e_above_mp_hull

assert df.e_above_mp_hull.isna().sum() == 0

e_above_mp_hull = df.e_above_mp_hull
target_col = "e_form_target"


mean = df.filter(like="pred").mean(axis=1) - df[target_col] + e_above_mp_hull

residual = mean - e_above_mp_hull

sort = np.argsort(e_above_mp_hull)

e_above_mp_hull = e_above_mp_hull[sort]
residual = residual[sort]

half_window = 0.02
increment = 0.002
bottom, top = -0.2, 0.3
bins = np.arange(bottom, top, increment)

rolling_means = np.zeros_like(bins)
rolling_stds = np.zeros_like(bins)

for j, b in enumerate(bins):
    low = b - half_window
    high = b + half_window

    mask = (e_above_mp_hull <= high) & (e_above_mp_hull > low)
    rolling_means[j] = residual.iloc[mask.values].abs().mean()
    rolling_stds[j] = std_err_of_mean(np.abs(residual.iloc[mask.values]))
fig, ax = plt.subplots(1, figsize=(10, 9))

ax.plot(bins, rolling_means)

ax.fill_between(
    bins, rolling_means + rolling_stds, rolling_means - rolling_stds, alpha=0.3
)

scale_bar = AnchoredSizeBar(
    ax.transData,
    2 * half_window,
    "40 meV",
    "lower left",
    pad=0,
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

arrowprops = dict(facecolor="black", width=0.5, headwidth=5, headlength=5)
ax.annotate(
    xy=(0.055, 0.05),
    xytext=(0.12, 0.05),
    arrowprops=arrowprops,
    text="Corrected\nGGA DFT\nAccuracy",
    verticalalignment="center",
    horizontalalignment="left",
)
ax.annotate(
    xy=(0.105, 0.1),
    xytext=(0.16, 0.1),
    arrowprops=arrowprops,
    text="GGA DFT\nAccuracy",
    verticalalignment="center",
    horizontalalignment="left",
)

ax.text(0, 0.13, r"$|\Delta E_{Hull-MP}| > $MAE", horizontalalignment="center")

ax.set(xlabel=r"$\Delta E_{Hull-MP}$ / eV per atom", ylabel="MAE / eV per atom")

ax.set(xlim=(bottom, top), ylim=(0.0, 0.14))

img_path = f"{PKG_DIR}/plots/{today}-moving-error-wbm-{rare=}.pdf"
# plt.savefig(img_path)

plt.show()
