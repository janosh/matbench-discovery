# %%
import pandas as pd

from matbench_discovery import ROOT, today
from matbench_discovery.plot_scripts import df_wbm
from matbench_discovery.plots import rolling_mae_vs_hull_dist

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"


# %%
data_path = (
    f"{ROOT}/data/2022-06-11-from-rhys/wren-mp-initial-structures.csv"
    # f"{ROOT}/models/wrenformer/2022-11-15-wrenformer-IS2RE-preds.csv"
)
df = pd.read_csv(data_path).set_index("material_id")
legend_label = "Wren"


# %%
df["e_above_hull_mp"] = df_wbm.e_above_hull_mp2020_corrected_ppd_mp

assert all(n_nans := df.isna().sum() == 0), f"Found {n_nans} NaNs"

target_col = "e_form_target"
# target_col = "e_form_per_atom"
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
    e_above_hull_true=df.e_above_hull_mp,
    label=legend_label,
)

fig = ax.figure
fig.set_size_inches(10, 9)
ax.legend(loc="lower right", frameon=False)

img_path = f"{ROOT}/figures/{today}-rolling-mae-vs-hull-dist.pdf"
# fig.savefig(img_path)
