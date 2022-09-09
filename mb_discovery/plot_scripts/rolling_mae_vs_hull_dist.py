# %%
from datetime import datetime

import pandas as pd

from mb_discovery import ROOT
from mb_discovery.plot_scripts import plt
from mb_discovery.plot_scripts.plot_funcs import rolling_mae_vs_hull_dist


__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

today = f"{datetime.now():%Y-%m-%d}"


# %%
markers = ["o", "v", "^", "H", "D", ""]

data_path = (
    f"{ROOT}/data/2022-06-11-from-rhys/wren-mp-initial-structures.csv"
    # f"{ROOT}/data/2022-08-16-wrenformer-preds.csv.bz2"
)
df = pd.read_csv(data_path).set_index("material_id")
legend_label = "Wren"
assert f"{legend_label.lower()}-" in data_path


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

assert all(n_nans := df.isna().sum() == 0), f"Found {n_nans} NaNs"

target_col = "e_form_target"
# --- or ---
# target_col = "e_form_per_atom_target"
# df["e_form_per_atom_target"] = df.e_form / df.n_sites

# make sure we average the expected number of ensemble member predictions
assert df.filter(regex=r"_pred_\d").shape[1] == 10

df["e_form_pres_ens"] = df.filter(regex=r"_pred_\d+").mean(axis=1)
df["e_above_hull_pred"] = df.e_form_pres_ens - df[target_col]


# %%
ax = rolling_mae_vs_hull_dist(
    e_above_hull_pred=df.e_above_hull_pred,
    e_above_hull_true=df.e_above_mp_hull,
    label=legend_label,
)

ax.figure.set_size_inches(10, 9)
ax.legend(loc="lower right", frameon=False)

img_path = f"{ROOT}/figures/{today}-rolling-mae-vs-hull-dist-{rare=}.pdf"
plt.savefig(img_path)
