"""Plot rolling MAE as a function of hull distance for a single model but split per WBM
batch in a single plot.
"""

# %%
import matplotlib.pyplot as plt
import pymatviz as pmv
from pymatviz.enums import Key
from pymatviz.utils import MATPLOTLIB, PLOTLY

from matbench_discovery import PDF_FIGS, SITE_FIGS, today
from matbench_discovery.data import Model
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.plots import rolling_mae_vs_hull_dist
from matbench_discovery.preds import df_each_pred, df_preds
from matbench_discovery.preds import models as all_models

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

batch_col = "batch_idx"
df_each_pred[batch_col] = "Batch " + df_each_pred.index.str.split("-").str[1]
df_err, df_std = None, None  # variables to cache rolling MAE and std
models = globals().get("models", all_models)


test_subset = globals().get("test_subset", TestSubset.uniq_protos)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(Key.uniq_proto)
    df_each_pred = df_each_pred.loc[df_preds.index]


# %% plotly version
for model in models:
    df_pivot = df_each_pred.pivot(columns=batch_col, values=model)  # noqa: PD010

    fig, df_err, df_std = rolling_mae_vs_hull_dist(
        e_above_hull_true=df_preds[MbdKey.each_true],
        e_above_hull_preds=df_pivot,
        # df_rolling_err=df_err,
        # df_err_std=df_std,
        backend=PLOTLY,
        show_dummy_mae=False,
        with_sem=False,
    )
    # if error is low, move legend to the top left
    leg_y = 1 if df_err.max().max() < 0.1 else 0.02
    y_anchor = "top" if leg_y == 1 else "bottom"
    fig.layout.legend.update(
        title=f"<b>{model}</b>",
        x=0.02,
        y=leg_y,
        bgcolor="rgba(0,0,0,0)",
        yanchor=y_anchor,
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
    pmv.save_fig(fig, f"{SITE_FIGS}/{img_path}.svelte")
    pmv.save_fig(fig, f"{PDF_FIGS}/{img_path}.pdf", width=500, height=330)


# %% matplotlib version
fig, ax = plt.subplots(1, figsize=(10, 9))
markers = ("o", "v", "^", "H", "D")
if len(markers) != 5:
    raise ValueError("Need 5 markers for 5 batches")
    # number of iterations of element substitution in WBM data set

model = Model.chgnet.label

for idx, marker in enumerate(markers, start=1):
    # select all rows from WBM step=idx
    df_step = df_preds[df_preds.index.str.startswith(f"wbm-{idx}-")]
    df_each_step = df_each_pred[df_each_pred.index.str.startswith(f"wbm-{idx}-")]

    title = f"Batch {idx} ({len(df_step.filter(like='e_').dropna()):,})"
    if not (1e4 < len(df_step) < 1e5):
        raise ValueError(f"WBM batches are 30k-50k in length, got {len(df_step)}")
    if any(df_step.index != df_each_step.index):
        raise ValueError("Index mismatch between df_step and df_each_step")

    ax, df_err, df_std = rolling_mae_vs_hull_dist(
        e_above_hull_true=df_step[MbdKey.each_true],
        e_above_hull_preds={title: df_each_step[model]},
        label=title,
        marker=marker,
        markevery=20,
        markerfacecolor="white",
        markeredgewidth=2.5,
        backend=MATPLOTLIB,  # don't change, code here not plotly compatible
        ax=ax,
        just_plot_lines=idx > 1,
        pbar=False,
    )


ax.legend(loc="lower right", frameon=False)
ax.set(title=f"{today} {model}")
for line in ax.lines:
    line._linewidth *= 3  # noqa: SLF001
    line.set_markersize(10)
