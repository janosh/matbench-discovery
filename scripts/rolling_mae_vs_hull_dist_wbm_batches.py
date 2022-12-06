# %%
from matbench_discovery import ROOT, today
from matbench_discovery.load_preds import load_df_wbm_with_preds
from matbench_discovery.plots import plt, rolling_mae_vs_hull_dist

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"


# %%
df_wbm = load_df_wbm_with_preds(models=["Wren", "Wrenformer"]).round(3)

e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"
target_col = "e_form_per_atom_mp2020_corrected"


# %%
model_name = "Wrenformer"
fig, ax = plt.subplots(1, figsize=(10, 9))
markers = ("o", "v", "^", "H", "D")
assert len(markers) == 5  # number of WBM rounds of element substitution

for idx, marker in enumerate(markers, 1):
    # select all rows from WBM step=idx
    df_step = df_wbm[df_wbm.index.str.startswith(f"wbm-step-{idx}")]

    title = f"Batch {idx} ({len(df_step.filter(like='e_').dropna()):,})"
    assert 1e4 < len(df_step) < 1e5, print(f"{len(df_step) = :,}")

    rolling_mae_vs_hull_dist(
        e_above_hull_true=df_step[e_above_hull_col],
        e_above_hull_error=df_step[target_col] - df_step[model_name],
        ax=ax,
        label=title,
        marker=marker,
        markevery=20,
        markerfacecolor="white",
        markeredgewidth=2.5,
    )


ax.legend(loc="lower right", frameon=False)
ax.set(title=f"{today} model={model_name}")


img_path = f"{ROOT}/figures/{today}-{model_name}-rolling-mae-vs-hull-dist-wbm-batches"
# fig.savefig(f"{img_path}.pdf")
