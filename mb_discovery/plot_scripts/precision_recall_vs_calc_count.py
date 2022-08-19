# %%
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from mb_discovery import ROOT


__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

today = f"{datetime.now():%Y-%m-%d}"

plt.rc("savefig", bbox="tight", dpi=200)
plt.rcParams["figure.constrained_layout.use"] = True
plt.rc("figure", dpi=150)
plt.rc("font", size=18)


# %%
df_hull = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")
dfs: dict[str, pd.DataFrame] = {}

for model_name in ("Wren", "CGCNN", "Voronoi"):
    dfs[model_name] = pd.read_csv(
        f"{ROOT}/data/2022-06-11-from-rhys/{model_name}-mp-initial-structures.csv"
    ).set_index("material_id")

dfs["m3gnet"] = pd.read_json(
    f"{ROOT}/data/2022-08-16-m3gnet-wbm-relax-results.json.gz"
).set_index("material_id")


# %%
fig, ax = plt.subplots(1, 1, figsize=(10, 9))

for model_name, color in zip(
    ("Wren", "CGCNN", "Voronoi", "M3GNet"),
    ("tab:blue", "tab:orange", "tab:red", "tab:green"),
):
    df = dfs[model_name]
    df = df.rename(columns={"e_form_wbm": "e_form_target"})

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

    if df.filter(regex=r"_pred_\d").shape[1] > 1:
        assert df.filter(regex=r"_pred_\d").shape[1] == 10

        model_preds = df.filter(regex=r"_pred_\d").mean(axis=1)

    elif model_name == "M3GNet":
        model_preds = df.e_form_m3gnet
    else:
        model_preds = df.e_form_pred

    residual = model_preds - df[target_col] + e_above_mp_hull

    # epistemic_var = df.filter(regex=r"_pred_\d").var(axis=1, ddof=0)

    # aleatoric_var = (df.filter(like="_ale_") ** 2).mean(axis=1)

    # full_std = (epistemic_var + aleatoric_var) ** 0.5

    # criterion = "std"
    # test = residual + full_std

    # criterion = "neg"
    # test = residual - full_std

    criterion = "energy"

    # thresh = 0.02
    thresh = 0
    # thresh = 0.10

    n_true_pos = len(
        e_above_mp_hull[(e_above_mp_hull <= thresh) & (residual <= thresh)]
    )
    n_false_neg = len(
        e_above_mp_hull[(e_above_mp_hull <= thresh) & (residual > thresh)]
    )

    n_total_pos = n_true_pos + n_false_neg

    sort = np.argsort(residual)
    e_above_mp_hull = e_above_mp_hull[sort]
    residual = residual[sort]

    e_type = "pred"
    true_pos_cumsum = ((e_above_mp_hull <= thresh) & (residual <= thresh)).cumsum()
    false_neg_cumsum = ((e_above_mp_hull <= thresh) & (residual > thresh)).cumsum()
    false_pos_cumsum = ((e_above_mp_hull > thresh) & (residual <= thresh)).cumsum()
    true_neg_cumsum = ((e_above_mp_hull > thresh) & (residual > thresh)).cumsum()
    xlabel = r"$\Delta E_{Hull-Pred}$ / eV per atom"

    ppv = true_pos_cumsum / (true_pos_cumsum + false_pos_cumsum) * 100
    tpr = true_pos_cumsum / n_total_pos * 100

    end = np.argmax(tpr)

    x = np.arange(len(ppv))[:end]

    precision_curve = interp1d(x, ppv[:end], kind="cubic")
    rolling_recall_curve = interp1d(x, tpr[:end], kind="cubic")

    line_kwargs = dict(
        linewidth=3, color=color, markevery=[-1], marker="x", markersize=14, mew=2.5
    )
    ax.plot(x[::100], precision_curve(x[::100]), linestyle="-", **line_kwargs)
    ax.plot(x[::100], rolling_recall_curve(x[::100]), linestyle=":", **line_kwargs)
    ax.plot((0, 0), (0, 0), label=model_name, **line_kwargs)


ax.set(xlabel="Number of Calculations", ylabel="Percentage")

ax.set(xlim=(0, 8e4), ylim=(0, 100))

model_legend = ax.legend(frameon=False, loc="lower right")
ax.add_artist(model_legend)

[precision] = ax.plot((0, 0), (0, 0), "black", linestyle="-")
[recall] = ax.plot((0, 0), (0, 0), "black", linestyle=":")
ax.legend(
    [precision, recall], ["Precision", "Recall"], frameon=False, loc="upper right"
)

img_path = (
    f"{ROOT}/figures/{today}-precision-recall-vs-calc-count-"
    f"{e_type=}-{criterion=}-{rare=}.pdf"
)
# plt.savefig(img_path)
