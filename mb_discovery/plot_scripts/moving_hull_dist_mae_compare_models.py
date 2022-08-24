# %%
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from scipy.stats import sem as std_err_of_mean

from mb_discovery import ROOT


__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

today = f"{datetime.now():%Y-%m-%d}"

plt.rc("savefig", bbox="tight", dpi=200)
plt.rcParams["figure.constrained_layout.use"] = True
plt.rc("figure", dpi=150)
plt.rc("font", size=16)


# %%
dfs: dict[str, pd.DataFrame] = {}
dfs["Wren"] = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wren-mp-initial-structures.csv"
).set_index("material_id")
dfs["CGCNN ISRE"] = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/cgcnn-mp-initial-structures.csv"
).set_index("material_id")
dfs["CGCNN Relaxed"] = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/cgcnn-mp-cse.csv"
).set_index("material_id")
dfs["Voronoi ISRE"] = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/voronoi-mp-initial-structures.csv"
).set_index("material_id")
dfs["Voronoi Relaxed"] = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/voronoi-mp-cse.csv"
).set_index("material_id")

df_hull = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")


# %%
fig, ax = plt.subplots(1, figsize=(10, 9))

half_window = 0.02
increment = 0.002
bottom, top = -0.2, 0.38
# bottom, top = -0.2, 0.3
target_col = "e_form_target"

for model_name, df in dfs.items():

    rare = "all"

    # from pymatgen.core import Composition
    # rare = "no-lanthanides"
    # df["contains_rare_earths"] = df.composition.map(
    #     lambda x: any(el.is_rare_earth_metal for el in Composition(x))
    # )
    # df = df.query("~contains_rare_earths")

    df["e_above_mp_hull"] = df_hull.e_above_mp_hull
    assert df.e_above_mp_hull.isna().sum() == 0

    # make sure we average the expected number of ensemble member predictions
    n_pred_cols = df.filter(like=r"_pred_").shape[1]
    if n_pred_cols > 1:
        assert n_pred_cols == 10, f"{n_pred_cols = }, expected 10"
        model_preds = df.filter(like=r"_pred_").mean(axis=1)
    else:
        model_preds = df.e_form_pred
    df["e_above_mp_hull_pred"] = model_preds - df[target_col] + df.e_above_mp_hull

    # epistemic_var = df.filter(like=r"_pred_").var(axis=1, ddof=0)
    # aleatoric_var = (df.filter(like="_ale_") ** 2).mean(axis=1)

    # std_total = (epistemic_var + aleatoric_var) ** 0.5

    df["residual"] = df.e_above_mp_hull_pred - df.e_above_mp_hull
    bins = np.arange(bottom, top, increment)

    rolling_maes = np.zeros_like(bins)
    rolling_mae_stds = np.zeros_like(bins)

    for idx, bin_center in enumerate(bins):
        low = bin_center - half_window
        high = bin_center + half_window

        mask = (df.e_above_mp_hull <= high) & (df.e_above_mp_hull > low)
        rolling_maes[idx] = df.residual[mask].abs().mean()
        rolling_mae_stds[idx] = std_err_of_mean(df.residual[mask].abs())

    print(model_name)
    print(f"  MAE: {df.residual.abs().mean():.3f}")
    print(f"  Min E: {rolling_maes.min():.3f}")
    assert rolling_maes.min() < df.residual.abs().mean() < rolling_maes.max()

    ax.plot(bins, rolling_maes, label=model_name)

    ax.fill_between(
        bins,
        rolling_maes + rolling_mae_stds,
        rolling_maes - rolling_mae_stds,
        alpha=0.3,
    )

scale_bar = AnchoredSizeBar(
    ax.transData,
    2 * half_window,
    "40 meV",
    "lower left",
    pad=0,
    borderpad=1,
    frameon=False,
    size_vertical=0.003,
)
ax.add_artist(scale_bar)

ax.plot((0.05, 0.5), (0.05, 0.5), color="grey", linestyle="--", alpha=0.3)
ax.plot((-0.5, -0.05), (0.5, 0.05), color="grey", linestyle="--", alpha=0.3)
ax.plot((-0.05, 0.05), (0.05, 0.05), color="grey", linestyle="--", alpha=0.3)
ax.plot((-0.1, 0.1), (0.1, 0.1), color="grey", linestyle="--", alpha=0.3)

ax.plot((0, 0.05), (0, 0.05), color="grey", linestyle="--", alpha=0.3)
ax.plot((-0.05, 0), (0.05, 0), color="grey", linestyle="--", alpha=0.3)

ax.set(xlabel=r"$\Delta E_{Hull-MP}$ / eV per atom", ylabel="MAE / eV per atom")
ax.set(xlim=(bottom, top), ylim=(0.0, 0.14))

ax.legend(frameon=False, loc="lower right")

img_path = f"{ROOT}/figures/{today}-moving-hull-dist-mae-compare-models-{rare=}.pdf"
# plt.savefig(img_path)

plt.show()
