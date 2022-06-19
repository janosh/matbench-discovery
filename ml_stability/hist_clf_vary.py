# %%
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from ml_stability import ROOT


__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

today = f"{datetime.now():%Y-%m-%d}"

plt.rc("font", size=18)
plt.rc("savefig", bbox="tight", dpi=200)
plt.rcParams["figure.constrained_layout.use"] = True
plt.rc("figure", dpi=150, titlesize=20)


# %%
fig, ax = plt.subplots(1, 1, figsize=(10, 9))

df_hull = pd.read_csv(f"{ROOT}/data/wbm_e_above_mp_hull.csv").set_index("material_id")


for model_name, color in zip(
    # ["wren", "cgcnn", "cgcnn-d"],
    # ["tab:blue", "tab:red", "tab:purple"],
    ["wren", "voronoi", "cgcnn"],
    ["tab:blue", "tab:orange", "tab:red"],
):
    df = pd.read_csv(f"{ROOT}/data/{model_name}-mp-initial-structures.csv").set_index(
        "material_id"
    )

    df["e_above_hull"] = df_hull.e_above_hull

    df = df.dropna(subset=["e_above_hull"])

    rare = "all"

    # rare = "nla"
    # df = df[
    #     ~df["composition"].apply(
    #         lambda x: any(el.is_rare_earth_metal for el in Composition(x).elements)
    #     )
    # ]

    e_above_hull = df.e_above_hull.to_numpy().ravel()

    # tar = df[tar_cols].to_numpy().ravel() - e_hull
    tar_f = df.filter(like="target").to_numpy().ravel()

    # mean = np.average(pred, axis=0) - e_hull
    mean = df.filter(like="pred").T.mean(axis=0) - tar_f + e_above_hull

    # epistemic_std = np.var(pred, axis=0, ddof=0)

    # aleatoric_std = np.mean(np.square(df.filter(like="ale")), axis=0)

    # full_std = np.sqrt(epistemic_std + aleatoric_std)

    # crit = "std"
    # test = mean + full_std

    # crit = "neg"
    # test = mean - full_std

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

    tp = len(e_above_hull[(e_above_hull <= thresh) & (test <= thresh)])
    fn = len(e_above_hull[(e_above_hull <= thresh) & (test > thresh)])

    pos = tp + fn

    sort = np.argsort(test)
    e_above_hull = e_above_hull[sort]
    test = test[sort]

    e_type = "pred"
    tp = np.asarray((e_above_hull <= thresh) & (test <= thresh))
    fn = np.asarray((e_above_hull <= thresh) & (test > thresh))
    fp = np.asarray((e_above_hull > thresh) & (test <= thresh))
    tn = np.asarray((e_above_hull > thresh) & (test > thresh))
    xlabel = r"$\Delta E_{Hull-Pred}$ / eV per atom"

    c_tp = np.cumsum(tp)
    c_fn = np.cumsum(fn)
    c_fp = np.cumsum(fp)
    c_tn = np.cumsum(tn)

    ppv = c_tp / (c_tp + c_fp) * 100
    tpr = c_tp / pos * 100

    end = np.argmax(tpr)

    x = np.arange(len(ppv))[:end]

    f_ppv = interp1d(x, ppv[:end], kind="cubic")
    f_tpr = interp1d(x, tpr[:end], kind="cubic")

    line_kwargs = dict(
        linewidth=2, color=color, markevery=[-1], marker="x", markersize=14, mew=2.5
    )
    ax.plot(x[::100], f_tpr(x[::100]), linestyle=":", **line_kwargs)

    ax.plot(x[::100], f_ppv(x[::100]), linestyle="-", **line_kwargs)


# ax.set_xticks((0, 2.5e4, 5e4, 7.5e4))
ax.set_xticks((0, 2e4, 4e4, 6e4, 8e4))

ax.set_ylabel("Percentage")
ax.set_xlabel("Number of Calculations")

ax.set_xlim((0, 8e4))
# ax.set_xlim((0, 75000))
ax.set_ylim((0, 100))

ax.plot((-1, -1), (-1, -1), color="tab:blue")
ax.plot((-1, -1), (-1, -1), color="tab:red")
ax.plot((-1, -1), (-1, -1), color="tab:orange")

# ax.plot((-1, -1), (-1, -1), color="tab:purple")

ax.plot((-1, -1), (-1, -1), "k", linestyle="-")
ax.plot((-1, -1), (-1, -1), "k", linestyle=":")

lines = ax.get_lines()

legend1 = ax.legend(
    lines[-2:], ["Precision", "Recall"], frameon=False, loc="upper right"
)
legend2 = ax.legend(
    lines[-5:-2],
    # ["Wren (This Work)", "CGCNN Pre-relax", "CGCNN-D Pre-relax"],
    ["Wren (This Work)", "CGCNN Pre-relax", "Voronoi Pre-relax"],
    frameon=False,
    loc="lower right",
)

ax.add_artist(legend1)
# plt.gca().add_artist(legend1)

ax.set_aspect(1.0 / ax.get_data_ratio())


# plt.savefig(f"{PKG_DIR}/plots/{today}-vary-{e_type}-{crit}-{rare}.pdf")
# # plt.savefig(f"{PKG_DIR}/plots/{today}-vary-{e_type}-{crit}-{rare}.png")

# plt.show()
