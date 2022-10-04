# %%
from datetime import datetime

import pandas as pd

from mb_discovery import ROOT
from mb_discovery.plots import plt, rolling_mae_vs_hull_dist

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

today = f"{datetime.now():%Y-%m-%d}"


# %%
rare = "all"

df_wren = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wren-mp-initial-structures.csv"
).set_index("material_id")

df_hull = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")

df_wren["e_above_mp_hull"] = df_hull.e_above_mp_hull
assert df_wren.e_above_mp_hull.isna().sum() == 0

target_col = "e_form_target"

# make sure we average the expected number of ensemble member predictions
assert df_wren.filter(regex=r"_pred_\d").shape[1] == 10

df_wren["e_above_hull_pred"] = (
    df_wren.filter(regex=r"_pred_\d").mean(axis=1) - df_wren[target_col]
)


# %%
fig, ax = plt.subplots(1, figsize=(10, 9))
markers = ("o", "v", "^", "H", "D")
assert len(markers) == 5  # number of WBM rounds of element substitution

for idx, marker in enumerate(markers, 1):
    df = df_wren[df_wren.index.str.startswith(f"wbm-step-{idx}")]
    title = f"Batch {idx} ({len(df.filter(like='e_').dropna()):,})"
    assert 1e4 < len(df) < 1e5, print(f"{len(df) = :,}")

    rolling_mae_vs_hull_dist(
        e_above_hull_pred=df.e_above_hull_pred,
        e_above_hull_true=df.e_above_mp_hull,
        ax=ax,
        label=title,
        marker=marker,
        markevery=20,
        markerfacecolor="white",
        markeredgewidth=2.5,
    )


ax.legend(loc="lower right", frameon=False)


img_path = f"{ROOT}/figures/{today}-rolling-mae-vs-hull-dist-wbm-batches-{rare=}.pdf"
# plt.savefig(img_path)
