# %%
import pandas as pd

from matbench_discovery import ROOT, today
from matbench_discovery.load_preds import df_wbm
from matbench_discovery.plots import plt, rolling_mae_vs_hull_dist

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"


# %%
df_wren = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wren-mp-initial-structures.csv"
).set_index("material_id")

df_wrenformer = pd.read_csv(
    f"{ROOT}/models/wrenformer/2022-11-15-wrenformer-IS2RE-preds.csv"
).set_index("material_id")


# %%
model_name = "wren"
df = {"wren": df_wren, "wrenformer": df_wrenformer}[model_name]

df["e_above_hull_mp"] = df_wbm.e_above_hull_mp2020_corrected_ppd_mp
assert df.e_above_hull_mp.isna().sum() == 0

possible_targets = (
    "e_form_per_atom_mp2020_corrected e_form_per_atom e_form_target".split()
)
target_col = next(filter(lambda x: x in df, possible_targets))

# make sure we average the expected number of ensemble member predictions
assert df.filter(regex=r"_pred_\d").shape[1] == 10

df["e_above_hull_pred"] = df.filter(regex=r"_pred_\d").mean(axis=1) - df[target_col]


# %%
fig, ax = plt.subplots(1, figsize=(10, 9))
markers = ("o", "v", "^", "H", "D")
assert len(markers) == 5  # number of WBM rounds of element substitution

for idx, marker in enumerate(markers, 1):
    # select all rows from WBM step=idx
    df_step = df[df.index.str.startswith(f"wbm-step-{idx}")]

    title = f"Batch {idx} ({len(df_step.filter(like='e_').dropna()):,})"
    assert 1e4 < len(df_step) < 1e5, print(f"{len(df_step) = :,}")

    rolling_mae_vs_hull_dist(
        e_above_hull_pred=df_step.e_above_hull_pred,
        e_above_hull_true=df_step.e_above_hull_mp,
        ax=ax,
        label=title,
        marker=marker,
        markevery=20,
        markerfacecolor="white",
        markeredgewidth=2.5,
    )


ax.legend(loc="lower right", frameon=False)
ax.set(title=f"{today} model={model_name}")


img_name = f"{today}-{model_name}-rolling-mae-vs-hull-dist-wbm-batches"
# fig.savefig(f"{ROOT}/figures/{img_name}.pdf")
