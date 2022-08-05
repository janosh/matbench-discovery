# %%
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from ml_stability import PKG_DIR, ROOT


__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

today = f"{datetime.now():%Y-%m-%d}"

plt.rc("font", size=18)
plt.rc("savefig", bbox="tight", dpi=200)
plt.rcParams["figure.constrained_layout.use"] = True
plt.rc("figure", dpi=150, titlesize=20)


# %%
df_hull = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")


# %%
fig, ax = plt.subplots(1, 1, figsize=(10, 9))

for model_name, color in zip(
    # ["wren", "cgcnn", "cgcnn-d"],
    # ["tab:blue", "tab:red", "tab:purple"],
    ["wren", "voronoi", "cgcnn"],
    ["tab:blue", "tab:orange", "tab:red"],
):
    data_path = (
        f"{ROOT}/data/2022-06-11-from-rhys/{model_name}-mp-initial-structures.csv"
    )
    df = pd.read_csv(data_path).set_index("material_id")

    df["e_above_mp_hull"] = df_hull.e_above_mp_hull

    assert df.e_above_mp_hull.isna().sum() == 0

    target_col = "e_form_target"
    rare = "all"

    # from pymatgen.core import Composition
    # rare = "no-lanthanides"
    # df["contains_rare_earths"] = df.composition.map(
    #     lambda x: any(el.is_rare_earth_metal for el in Composition(x))
    # )
    # df = df.query("~contains_rare_earths")

    e_above_mp_hull = df.e_above_mp_hull

    # mean = df.filter(like="pred").mean(axis=1) - e_hull
    mean = df.filter(like="pred").mean(axis=1) - df[target_col] + e_above_mp_hull

    # epistemic_var = df.filter(like="pred").var(axis=1, ddof=0)

    # aleatoric_var = (df.filter(like="ale") ** 2).mean(axis=1)

    # full_std = (epistemic_var + aleatoric_var) ** 0.5

    # criterion = "std"
    # test = mean + full_std

    # criterion = "neg"
    # test = mean - full_std

    criterion = "energy"

    bins = 200
    # xlim = (-0.2, 0.2)
    xlim = (-0.4, 0.4)
    # xlim = (-1, 1)

    # thresh = 0.02
    thresh = 0.00
    # thresh = 0.10
    xticks = (-0.4, -0.2, 0, 0.2, 0.4)
    # yticks = (0, 300, 600, 900, 1200)

    n_true_pos = len(e_above_mp_hull[(e_above_mp_hull <= thresh) & (mean <= thresh)])
    n_false_neg = len(e_above_mp_hull[(e_above_mp_hull <= thresh) & (mean > thresh)])

    n_total_pos = n_true_pos + n_false_neg

    sort = np.argsort(mean)
    e_above_mp_hull = e_above_mp_hull[sort]
    mean = mean[sort]

    e_type = "pred"
    true_pos_cumsum = ((e_above_mp_hull <= thresh) & (mean <= thresh)).cumsum()
    false_neg_cumsum = ((e_above_mp_hull <= thresh) & (mean > thresh)).cumsum()
    false_pos_cumsum = ((e_above_mp_hull > thresh) & (mean <= thresh)).cumsum()
    true_neg_cumsum = ((e_above_mp_hull > thresh) & (mean > thresh)).cumsum()
    xlabel = r"$\Delta E_{Hull-Pred}$ / eV per atom"

    ppv = true_pos_cumsum / (true_pos_cumsum + false_pos_cumsum) * 100
    tpr = true_pos_cumsum / n_total_pos * 100

    end = np.argmax(tpr)

    x = np.arange(len(ppv))[:end]

    precision_curve = interp1d(x, ppv[:end], kind="cubic")
    rolling_recall_curve = interp1d(x, tpr[:end], kind="cubic")

    line_kwargs = dict(
        linewidth=2, color=color, markevery=[-1], marker="x", markersize=14, mew=2.5
    )
    ax.plot(x[::100], precision_curve(x[::100]), linestyle="-", **line_kwargs)
    ax.plot(x[::100], rolling_recall_curve(x[::100]), linestyle=":", **line_kwargs)


ax.set(xlabel="Number of Calculations", ylabel="Percentage")

xlim = (0, 8e4)
ax.set(xlim=xlim, ylim=(0, 100), xticks=np.linspace(*xlim, 5))

ax.plot((0, 0), (0, 0), color="tab:blue", label="Wren (This Work)")
ax.plot((0, 0), (0, 0), color="tab:red", label="CGCNN Pre-relax")
ax.plot((0, 0), (0, 0), color="tab:orange", label="Voronoi Pre-relax")
legend_1 = ax.legend(frameon=False, loc="lower right")
ax.add_artist(legend_1)

[prec] = ax.plot((0, 0), (0, 0), "black", linestyle="-")
[recall] = ax.plot((0, 0), (0, 0), "black", linestyle=":")
ax.legend([prec, recall], ["Precision", "Recall"], frameon=False, loc="upper right")

img_path = f"{PKG_DIR}/plots/{today}-vary-{e_type=}-{criterion=}-{rare=}.pdf"
# plt.savefig(img_path)
