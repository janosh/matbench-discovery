# %%
from datetime import datetime

import matplotlib.pyplot as plt
import pandas as pd

from mb_discovery import ROOT
from mb_discovery.plot_scripts.plot_funcs import rolling_mae_vs_hull_dist


__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

today = f"{datetime.now():%Y-%m-%d}"

plt.rc("savefig", bbox="tight", dpi=200)
plt.rcParams["figure.constrained_layout.use"] = True
plt.rc("figure", dpi=200)
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
    assert df.isna().sum().sum() == 0

    # make sure we average the expected number of ensemble member predictions
    n_pred_cols = df.filter(like=r"_pred_").shape[1]
    if n_pred_cols > 1:
        assert n_pred_cols == 10, f"{n_pred_cols = }, expected 10"
        model_preds = df.filter(like=r"_pred_").mean(axis=1)
    else:
        model_preds = df.e_form_pred
    df["e_above_mp_hull_pred"] = model_preds - df[target_col] + df.e_above_mp_hull
    df["residual"] = df.e_above_mp_hull_pred - df.e_above_mp_hull

    rolling_mae_vs_hull_dist(
        df,
        e_above_hull_col="e_above_mp_hull",
        residual_col="residual",
        ax=ax,
        label=model_name,
    )

ax.legend(frameon=False, loc="lower right")

img_path = f"{ROOT}/figures/{today}-rolling-mae-vs-hull-dist-compare-models-{rare=}.pdf"
# plt.savefig(img_path)

plt.show()
