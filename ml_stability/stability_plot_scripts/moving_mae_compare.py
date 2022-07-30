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
df_cgcnn_pre = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/cgcnn-mp-initial-structures.csv"
).set_index("material_id")
df_cgcnn_rel = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/cgcnn-mp-cse.csv",
    # f"aviary/results/manuscript/wbm-step-{i+offsets}_cgcnn-pre_org.csv"
).set_index("material_id")
# df_cgcnn_dis = pd.read_csv(
#     f"aviary/results/manuscript/wbm-step-{i+offsets}_cgcnn-d_org.csv",
#     # f"aviary/results/manuscript/wbm-step-{i+offsets}_cgcnn-pre_org.csv",
# ).set_index("material_id")
df_vt_pre = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/voronoi-mp-initial-structures.csv"
).set_index("material_id")
df_vt_rel = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/voronoi-mp-cse.csv"
).set_index("material_id")
df_wren = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wren-mp-initial-structures.csv"
).set_index("material_id")


# %% Find MAD Voronoi
# mat_ids = list(
#     set(df_vt_pre.material_id.to_list()).intersection(df_vt_rel.material_id.to_list())
# )

# df_vt_pre = df_vt_pre.set_index("material_id")
# df_vt_rel = df_vt_rel.set_index("material_id")

# print(np.abs(
#         np.average(df_vt_pre.filter(like="pred").loc[mat_ids].to_numpy().T, axis=0) -
#         np.average(df_vt_rel.filter(like="pred").loc[mat_ids].to_numpy().T, axis=0)
#     ).mean()
# )
# print(
#     (
#         np.average(df_vt_pre.filter(like="pred").loc[mat_ids].to_numpy().T, axis=0) -
#         np.average(df_vt_rel.filter(like="pred").loc[mat_ids].to_numpy().T, axis=0)
#     ).mean()
# )


# %%
df_hull = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")

fig, ax = plt.subplots(1, figsize=(10, 9))

for df, model_name, linestyle in zip(
    # (df_wren, df_vt_pre, df_vt_rel, df_cgcnn_pre, df_cgcnn_dis, df_cgcnn_rel),
    # (
    #     "Wren (This Work)",
    #     "Voronoi Pre-relax",
    #     "Voronoi Relaxed",
    #     "CGCNN Pre-relax",
    #     "CGCNN-D Pre-relax",
    #     "CGCNN Relaxed",
    # ),
    (df_wren, df_vt_pre, df_vt_rel, df_cgcnn_pre, df_cgcnn_rel),
    (
        "Wren (This Work)",
        "Voronoi Pre-relax",
        "Voronoi Relaxed",
        "CGCNN Pre-relax",
        "CGCNN Relaxed",
    ),
    # (df_wren, df_vt_pre, df_vt_rel),
    # ("Wren (This Work)", "Voronoi Pre-relax", "Voronoi Relaxed"),
    # (df_wren, df_cgcnn_pre, df_cgcnn_rel),
    # ("Wren (This Work)", "CGCNN Pre-relax", "CGCNN Relaxed"),
    ("-", "--", ":", "-.", ":", "--"),
):

    rare = "all"

    # from pymatgen.core import Composition
    # rare = "no-lanthanides"
    # df["contains_rare_earths"] = df.composition.map(
    #     lambda x: any(el.is_rare_earth_metal for el in Composition(x))
    # )
    # df = df.query("~contains_rare_earths")

    df["e_above_hull"] = df_hull.e_above_hull
    df = df.dropna(subset=["e_above_hull"])
    tar = df.e_above_hull.to_numpy().ravel()

    tar_f = df.filter(like="target").to_numpy().ravel()

    pred = df.filter(like="pred").to_numpy().T
    mean = np.average(pred, axis=0) - tar_f + tar

    # epistemic_std = np.var(pred, axis=0, ddof=0)
    # aleatoric_std = np.mean(np.square(df.filter(like="ale")), axis=0)

    # full_std = np.sqrt(epistemic_std + aleatoric_std)

    res = mean - tar

    half_window = 0.02
    increment = 0.002
    bot, top = -0.2, 0.38
    # bot, top = -0.2, 0.3
    bins = np.arange(bot, top, increment)

    means = np.zeros_like(bins)
    std = np.zeros_like(bins)

    for j, b in enumerate(bins):
        low = b - half_window
        high = b + half_window

        means[j] = np.mean(np.abs(res[np.argwhere((tar <= high) & (tar > low))]))
        std[j] = std_err_of_mean(np.abs(res[np.argwhere((tar <= high) & (tar > low))]))

    print(model_name)
    print(f"  MAE: {abs(res).mean():.4f}")
    print(f"  Min E: {means.min():.4f}")

    ax.plot(bins, means, linestyle=linestyle, label=model_name)

    ax.fill_between(bins, means + std, means - std, alpha=0.3)

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
)

ax.add_artist(scalebar)

ax.plot((0.05, 0.5), (0.05, 0.5), color="grey", linestyle="--", alpha=0.3)
ax.plot((-0.5, -0.05), (0.5, 0.05), color="grey", linestyle="--", alpha=0.3)
ax.plot((-0.05, 0.05), (0.05, 0.05), color="grey", linestyle="--", alpha=0.3)
ax.plot((-0.1, 0.1), (0.1, 0.1), color="grey", linestyle="--", alpha=0.3)

ax.plot((0, 0.05), (0, 0.05), color="grey", linestyle="--", alpha=0.3)
ax.plot((-0.05, 0), (0.05, 0), color="grey", linestyle="--", alpha=0.3)

ax.set(xlabel=r"$\Delta E_{Hull-MP}$ / eV per atom", ylabel="MAE / eV per atom")
ax.set(xlim=(bot, top), ylim=(0.0, 0.14))

ax.legend(
    frameon=False,
    loc="lower right",
    # facecolor="white",
    # framealpha=1.0,
    # edgecolor="white",
)


plt.savefig(f"{PKG_DIR}/plots/{today}-moving-error-wbm-{rare}-compare.pdf")
# plt.savefig(f"{PKG_DIR}/plots/{today}-moving-error-wbm-{rare}-compare.png")
# plt.savefig(f"examples/plots/pdf/moving-error-wbm-{rare}-all.png")

plt.show()
