# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from scipy.stats import sem


plt.rcParams.update({"font.size": 20})

plt.rcParams["axes.linewidth"] = 2.5
plt.rcParams["lines.linewidth"] = 3.5
plt.rcParams["xtick.major.size"] = 7
plt.rcParams["xtick.major.width"] = 2.5
plt.rcParams["xtick.minor.size"] = 5
plt.rcParams["xtick.minor.width"] = 2.5
plt.rcParams["ytick.major.size"] = 7
plt.rcParams["ytick.major.width"] = 2.5
plt.rcParams["ytick.minor.size"] = 5
plt.rcParams["ytick.minor.width"] = 2.5
plt.rcParams["legend.fontsize"] = 20


# %%
fig, ax = plt.subplots(1, figsize=(10, 9))

markers = [
    "o",
    "v",
    "^",
    "H",
    "D",
    "",
]

df = pd.read_csv(
    f"/home/reag2/PhD/aviary/examples/manuscript/new_figs/wren-mp-init.csv",
    comment="#",
    na_filter=False,
)


rare = "all"
# rare = "nla"
# df = df[
#     ~df["composition"].apply(
#         lambda x: any(el.is_rare_earth_metal for el in Composition(x).elements)
#     )
# ]

df_hull = pd.read_csv(
    f"/home/reag2/PhD/aviary/examples/manuscript/new_figs/wbm_e_above_mp.csv",
    comment="#",
    na_filter=False,
)

df["E_hull"] = pd.to_numeric(
    df["material_id"].map(dict(zip(df_hull.material_id, df_hull.E_above_hull)))
)

df = df.dropna(axis=0, subset=["E_hull"])

tar = df["E_hull"].to_numpy().ravel()

tar_cols = [col for col in df.columns if "target" in col]
tar_f = df[tar_cols].to_numpy().ravel()

pred_cols = [col for col in df.columns if "pred" in col]
pred = df[pred_cols].to_numpy().T
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
    std[j] = sem(np.abs(res[np.argwhere((tar <= high) & (tar > low))]))

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

ineq = "|" + r"$\Delta$" + r"$\it{E}$" + r"$_{Hull-MP}$" + "| > MAE"

ax.text(0, 0.13, ineq, horizontalalignment="center")

ax.set_ylabel("MAE / eV per atom")
x_lab = r"$\Delta$" + r"$\it{E}$" + r"$_{Hull-MP}$" + " / eV per atom"

ax.set_xlabel(x_lab)

ax.set_ylim((0.0, 0.14))
ax.set_xlim((bot, top))

ax.set_aspect(1.0 / ax.get_data_ratio())

plt.tight_layout()

plt.savefig(f"examples/manuscript/new_figs/moving-error-wbm-{rare}.pdf")

plt.show()
