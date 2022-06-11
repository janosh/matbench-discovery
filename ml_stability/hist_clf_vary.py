# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


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

fig, ax = plt.subplots(1, 1, figsize=(10, 9))

df_hull = pd.read_csv(
    f"/home/reag2/PhD/aviary/examples/manuscript/new_figs/wbm_e_above_mp.csv",
    comment="#",
    na_filter=False,
)

e_hull_dict = dict(zip(df_hull.material_id, df_hull.E_above_hull))

for name, c, a in zip(
    # ["wren", "cgcnn", "cgcnn-d"],
    # ["tab:blue", "tab:red", "tab:purple"],
    ["wren", "voro", "cgcnn"],
    ["tab:blue", "tab:orange", "tab:red"],
    [1, 0.8, 0.8],
    # ["wren", "cgcnn"],
    # ["tab:blue", "tab:red"],
    # [1, 0.8],
):
    df = pd.read_csv(
        f"/home/reag2/PhD/aviary/examples/manuscript/new_figs/{name}-mp-init.csv",
        comment="#",
        na_filter=False,
    )

    df["E_hull"] = pd.to_numeric(df["material_id"].map(e_hull_dict))

    df = df.dropna(axis=0, subset=["E_hull"])

    init = len(df)

    rare = "all"

    # rare = "nla"
    # df = df[
    #     ~df["composition"].apply(
    #         lambda x: any(el.is_rare_earth_metal for el in Composition(x).elements)
    #     )
    # ]

    # print(1-len(df)/init)

    tar = df["E_hull"].to_numpy().ravel()

    print(len(tar))

    tar_cols = [col for col in df.columns if "target" in col]
    # tar = df[tar_cols].to_numpy().ravel() - e_hull
    tar_f = df[tar_cols].to_numpy().ravel()

    pred_cols = [col for col in df.columns if "pred" in col]
    pred = df[pred_cols].to_numpy().T
    # mean = np.average(pred, axis=0) - e_hull
    mean = np.average(pred, axis=0) - tar_f + tar

    epi = np.var(pred, axis=0, ddof=0)

    ale_cols = [col for col in df.columns if "ale" in col]
    if len(ale_cols) > 0:
        ales = df[ale_cols].to_numpy().T
        ale = np.mean(np.square(ales), axis=0)
    else:
        ale = 0

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

    alpha = 0.5
    # thresh = 0.02
    thresh = 0.00
    # thresh = 0.10
    xticks = (-0.4, -0.2, 0, 0.2, 0.4)
    # yticks = (0, 300, 600, 900, 1200)

    tp = len(tar[(tar <= thresh) & (test <= thresh)])
    fn = len(tar[(tar <= thresh) & (test > thresh)])

    pos = tp + fn
    null = pos / len(tar)

    sort = np.argsort(test)
    tar = tar[sort]
    test = test[sort]

    e_type = "pred"
    tp = np.asarray((tar <= thresh) & (test <= thresh))
    fn = np.asarray((tar <= thresh) & (test > thresh))
    fp = np.asarray((tar > thresh) & (test <= thresh))
    tn = np.asarray((tar > thresh) & (test > thresh))
    xlabel = (
        r"$\Delta$" + r"$\it{E}$" + r"$_{Hull-Pred}$" + " / eV per atom"
    )  # r"$\/(\frac{eV}{atom})$"

    # %%

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

    ax.plot(
        x[::100],
        f_tpr(x[::100]),
        linestyle=":",
        color=c,
        alpha=a,
        markevery=[-1],
        marker="x",
        markersize=14,
        mew=2.5,
    )

    ax.plot(
        x[::100],
        f_ppv(x[::100]),
        linestyle="-",
        color=c,
        alpha=a,
        markevery=[-1],
        marker="x",
        markersize=14,
        mew=2.5,
    )


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


fig.tight_layout()
plt.savefig(f"examples/manuscript/new_figs/vary-{e_type}-{crit}-{rare}.pdf")
# plt.savefig(f"examples/manuscript/pdf/vary-{e_type}-{crit}-{rare}.png")

plt.show()
