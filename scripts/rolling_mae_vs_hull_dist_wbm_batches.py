# %%
from matbench_discovery import FIGS, today
from matbench_discovery.plots import plt, rolling_mae_vs_hull_dist
from matbench_discovery.preds import df_wbm, e_form_col, each_true_col

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"


# %%
model_name = "Wrenformer"
fig, ax = plt.subplots(1, figsize=(10, 9))
markers = ("o", "v", "^", "H", "D")
assert len(markers) == 5  # number of iterations of element substitution in WBM data set

for idx, marker in enumerate(markers, 1):
    # select all rows from WBM step=idx
    df_step = df_wbm[df_wbm.index.str.startswith(f"wbm-{idx}-")]

    title = f"Batch {idx} ({len(df_step.filter(like='e_').dropna()):,})"
    assert 1e4 < len(df_step) < 1e5, print(f"{len(df_step) = :,}")

    ax, df_err, df_std = rolling_mae_vs_hull_dist(
        e_above_hull_true=df_step[each_true_col],
        e_above_hull_errors={title: df_step[e_form_col] - df_step[model_name]},
        label=title,
        marker=marker,
        markevery=20,
        markerfacecolor="white",
        markeredgewidth=2.5,
        backend="matplotlib",
        ax=ax,
        just_plot_lines=idx > 1,
    )


ax.legend(loc="lower right", frameon=False)
ax.set(title=f"{today} {model_name}")
for line in ax.lines:
    line._linewidth *= 3
    line.set_markersize(10)


img_path = f"{FIGS}/{today}-{model_name}-rolling-mae-vs-hull-dist-wbm-batches"
# fig.savefig(f"{img_path}.pdf")
