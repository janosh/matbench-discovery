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
markers = [
    "o",
    "v",
    "^",
    "H",
    "D",
    "",
]

tars = []
df_hull_list = []
df_list_wren = []
df_list_cgcnn_pre = []
df_list_cgcnn_rel = []
df_list_cgcnn_dis = []
df_list_vt_pre = []
df_list_vt_rel = []


df_cgcnn_pre = pd.read_csv(
    "/home/reag2/PhD/aviary/examples/manuscript/new_figs/cgcnn-mp-init.csv",
    comment="#",
    na_filter=False,
    # index_col="material_id"
)
df_cgcnn_rel = pd.read_csv(
    "/home/reag2/PhD/aviary/examples/manuscript/new_figs/cgcnn-mp-cse.csv",
    # f"/home/reag2/PhD/aviary/results/manuscript/step_{i+offsets}_cgcnn-pre_org.csv",
    comment="#",
    na_filter=False,
    # index_col="material_id"
)
# df_cgcnn_dis = pd.read_csv(
#     f"/home/reag2/PhD/aviary/results/manuscript/step_{i+offsets}_cgcnn-d_org.csv",
#     # f"/home/reag2/PhD/aviary/results/manuscript/step_{i+offsets}_cgcnn-pre_org.csv",
#     comment="#",
#     na_filter=False,
#     # index_col="material_id"
# )
df_vt_pre = pd.read_csv(
    "/home/reag2/PhD/aviary/examples/manuscript/new_figs/voro-mp-init.csv",
    comment="#",
    na_filter=False,
    # index_col="material_id"
)
df_vt_rel = pd.read_csv(
    "/home/reag2/PhD/aviary/examples/manuscript/new_figs/voro-mp-cse.csv",
    comment="#",
    na_filter=False,
    # index_col="material_id"
)
df_wren = pd.read_csv(
    "/home/reag2/PhD/aviary/examples/manuscript/new_figs/wren-mp-init.csv",
    comment="#",
    na_filter=False,
)


# %%
# Find MAD voro

# pred_cols = [col for col in df_vt_pre.columns if "pred" in col]

# mat_ids = list(set(df_vt_pre.material_id.to_list()).intersection(df_vt_rel.material_id.to_list()))

# df_vt_pre = df_vt_pre.set_index("material_id")
# df_vt_rel = df_vt_rel.set_index("material_id")

# print(np.abs(
#         np.average(df_vt_pre[pred_cols].loc[mat_ids].to_numpy().T, axis=0) -
#         np.average(df_vt_rel[pred_cols].loc[mat_ids].to_numpy().T, axis=0)
#     ).mean()
# )
# print(
#     (
#         np.average(df_vt_pre[pred_cols].loc[mat_ids].to_numpy().T, axis=0) -
#         np.average(df_vt_rel[pred_cols].loc[mat_ids].to_numpy().T, axis=0)
#     ).mean()
# )


# %%
df_hull = pd.read_csv(
    "/home/reag2/PhD/aviary/examples/manuscript/new_figs/wbm_e_above_mp.csv",
    comment="#",
    na_filter=False,
)

e_hull_dict = dict(zip(df_hull.material_id, df_hull.E_above_hull))

fig, ax = plt.subplots(1, figsize=(10, 9))

for df, n, l, a in zip(
    # (df_wren, df_vt_pre, df_vt_rel, df_cgcnn_pre, df_cgcnn_dis, df_cgcnn_rel),
    # ("Wren (This Work)", "Voronoi Pre-relax", "Voronoi Relaxed", "CGCNN Pre-relax", "CGCNN-D Pre-relax", "CGCNN Relaxed"),
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
    (1.0, 0.8, 0.8, 0.8, 0.8, 0.8),
):

    rare = "all"

    # rare = "nla"
    # df = df[
    #     ~df["composition"].apply(
    #         lambda x: any(el.is_rare_earth_metal for el in Composition(x).elements)
    #     )
    # ]

    df["E_hull"] = pd.to_numeric(df["material_id"].map(e_hull_dict))
    df = df.dropna(axis=0, subset=["E_hull"])
    tar = df["E_hull"].to_numpy().ravel()

    tar_cols = [col for col in df.columns if "target" in col]
    tar_f = df[tar_cols].to_numpy().ravel()

    pred_cols = [col for col in df.columns if "pred" in col]
    pred = df[pred_cols].to_numpy().T
    mean = np.average(pred, axis=0) - tar_f + tar

    epi = np.var(pred, axis=0, ddof=0)

    ale_cols = [col for col in df.columns if "ale" in col]
    if len(ale_cols) > 0:
        ales = df[ale_cols].to_numpy().T
        ale = np.mean(np.square(ales), axis=0)
    else:
        ale = 0

    both = np.sqrt(epi + ale)

    res = mean - tar

    print(f"{n}")
    print(f"MAE: {np.mean(np.abs(res))}")

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
        std[j] = sem(np.abs(res[np.argwhere((tar <= high) & (tar > low))]))

    print(f"Min E: {np.min(means)}")

    ax.plot(
        bins,
        means,
        linestyle=l,
        alpha=a,
        label=n,
    )

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

ax.plot((0, 0.05), (0, 0.05), color="grey", linestyle="--", alpha=0.3)
ax.plot((-0.05, 0), (0.05, 0), color="grey", linestyle="--", alpha=0.3)

ax.set_ylabel("MAE / eV per atom")
x_lab = r"$\Delta$" + r"$\it{E}$" + r"$_{Hull-MP}$" + " / eV per atom"

ax.set_xlabel(x_lab)

ax.set_ylim((0.0, 0.14))
ax.set_xlim((bot, top))
ax.legend(
    frameon=False,
    loc="lower right",
    # facecolor="white",
    # framealpha=1.0,
    # edgecolor="white",
)

ax.set_aspect(1.0 / ax.get_data_ratio())

plt.tight_layout()

plt.savefig(f"examples/manuscript/new_figs/moving-error-wbm-{rare}-compare.pdf")
# plt.savefig(f"examples/manuscript/pdf/moving-error-wbm-{rare}-compare.png")
# plt.savefig(f"examples/plots/pdf/moving-error-wbm-{rare}-all.png")

plt.show()
