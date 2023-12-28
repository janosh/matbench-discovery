"""Plot rolling MAE as a function of hull distance for a single model but split per WBM
batch in a single plot.
"""


# %%
from pymatviz.io import save_fig

from matbench_discovery import PDF_FIGS, SITE_FIGS, today
from matbench_discovery.plots import plt, rolling_mae_vs_hull_dist
from matbench_discovery.preds import df_each_pred, df_preds, e_form_col, each_true_col
from matbench_discovery.preds import models as all_models

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

batch_col = "batch_idx"
df_each_pred[batch_col] = "Batch " + df_each_pred.index.str.split("-").str[1]
df_err, df_std = None, None  # variables to cache rolling MAE and std
models = globals().get("models", all_models)


# %% matplotlib version
fig, ax = plt.subplots(1, figsize=(10, 9))
markers = ("o", "v", "^", "H", "D")
assert len(markers) == 5  # number of iterations of element substitution in WBM data set
model = "CHGNet"

for idx, marker in enumerate(markers, 1):
    # select all rows from WBM step=idx
    df_step = df_preds[df_preds.index.str.startswith(f"wbm-{idx}-")]

    title = f"Batch {idx} ({len(df_step.filter(like='e_').dropna()):,})"
    assert 1e4 < len(df_step) < 1e5, print(f"{len(df_step) = :,}")

    ax, df_err, df_std = rolling_mae_vs_hull_dist(
        e_above_hull_true=df_step[each_true_col],
        e_above_hull_errors={title: df_step[e_form_col] - df_step[model]},
        label=title,
        marker=marker,
        markevery=20,
        markerfacecolor="white",
        markeredgewidth=2.5,
        backend="matplotlib",  # don't change, code here not plotly compatible
        ax=ax,
        just_plot_lines=idx > 1,
        pbar=False,
    )


ax.legend(loc="lower right", frameon=False)
ax.set(title=f"{today} {model}")
for line in ax.lines:
    line._linewidth *= 3  # noqa: SLF001
    line.set_markersize(10)


# %% plotly version
for model in models:
    df_pivot = df_each_pred.pivot(columns=batch_col, values=model)

    fig, df_err, df_std = rolling_mae_vs_hull_dist(
        e_above_hull_true=df_preds[each_true_col],
        e_above_hull_errors=df_pivot,
        # df_rolling_err=df_err,
        # df_err_std=df_std,
        backend="plotly",
        show_dummy_mae=False,
        with_sem=False,
    )
    fig.layout.legend.update(
        title=f"<b>{model}</b>", x=0.02, y=0.02, bgcolor="rgba(0,0,0,0)"
    )
    fig.layout.margin.update(l=10, r=10, b=10, t=10)
    fig.layout.update(hovermode="x unified", hoverlabel_bgcolor="black")
    fig.update_traces(
        hovertemplate="y=%{y:.3f} eV",
        selector=lambda trace: trace.name.startswith("Batch"),
    )
    fig.show()

    model_snake_case = model.lower().replace(" + ", "-").replace(" ", "-")
    img_path = f"rolling-mae-vs-hull-dist-wbm-batches-{model_snake_case}"
    save_fig(fig, f"{SITE_FIGS}/{img_path}.svelte")
    save_fig(fig, f"{PDF_FIGS}/{img_path}.pdf", width=500, height=330)
