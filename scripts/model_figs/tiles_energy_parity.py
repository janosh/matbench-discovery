"""Parity plot of actual vs predicted e_above_hull and e_form_per_atom for all
models. Unlike energy_parity_traces.py, this script splits each model into a
separate subplot.
"""

# %%
import math

import numpy as np
import plotly.express as px
import pymatviz as pmv
from pymatviz.enums import Key
from pymatviz.utils import bin_df_cols

from matbench_discovery import PDF_FIGS, SITE_FIGS
from matbench_discovery.cli import cli_args
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey
from matbench_discovery.metrics.discovery import dfs_metrics
from matbench_discovery.plots import clf_colors

__author__ = "Janosh Riebesell"
__date__ = "2022-11-28"

if cli_args.energy_type == Key.each:
    e_pred_col = Key.each_pred
    e_true_col = MbdKey.each_true
elif cli_args.energy_type == Key.e_form:
    e_true_col = MbdKey.e_form_dft
    e_pred_col = Key.e_form_pred
else:
    raise ValueError(f"Unexpected {cli_args.energy_type=}")

# Get list of models from command line args, fall back to all models
models_to_plot = cli_args.models
test_subset = cli_args.test_subset

# Load predictions for specified models
df_preds = load_df_wbm_with_preds(
    models=models_to_plot, subset=cli_args.test_subset
).round(3)


# %%
facet_col = "Model"
hover_cols = (MbdKey.each_true, Key.formula)

df_melt = df_preds.melt(
    id_vars=(df_preds.index.name, MbdKey.e_form_dft, *hover_cols),
    var_name=facet_col,
    value_vars=[model.label for model in models_to_plot],
    value_name=Key.e_form_pred,
)

df_melt[Key.each_pred] = (
    df_melt[MbdKey.each_true] + df_melt[Key.e_form_pred] - df_melt[MbdKey.e_form_dft]
)

df_bin = bin_df_cols(
    df_melt,
    bin_by_cols=[e_true_col, e_pred_col],
    group_by_cols=[facet_col],
    n_bins=300,
    bin_counts_col=(bin_cnt_col := "bin counts"),
)
df_bin = df_bin.reset_index()

# sort legend and facet plots by MAE
legend_order = list(dfs_metrics[test_subset].T.MAE.sort_values().index)


# determine each point's classification to color them by
# now unused, can be used to color points by TP/FP/TN/FN
# true_pos, false_neg, false_pos, true_neg = classify_stable(
#     df_bin[e_true_col], df_bin[e_pred_col]
# )
# clf_col = "classified"
# df_bin[clf_col] = np.array(clf_labels)[
#     true_pos * 0 + false_neg * 1 + false_pos * 2 + true_neg * 3
# ]


# %% parity plot of DFT vs predicted hull distance with each model in separate subplot
models_to_plot = [
    model
    for model in models_to_plot
    if model.is_complete and (cli_args.show_non_compliant or model.is_compliant)
]

log_bin_cnt_col = f"log {bin_cnt_col}"
df_bin[log_bin_cnt_col] = np.log1p(df_bin[bin_cnt_col]).round(2)

n_cols = 3
if cli_args.use_full_rows:
    # drop last models that don't fit in last row
    n_rows = len(models_to_plot) // n_cols
    models_to_plot = models_to_plot[: n_rows * n_cols]
else:
    n_rows = math.ceil(len(models_to_plot) / n_cols)

fig = px.scatter(
    df_bin.query(f"{facet_col} in {[m.label for m in models_to_plot]}"),
    x=e_true_col,
    y=e_pred_col,
    facet_col=facet_col,
    facet_col_wrap=n_cols,
    color=log_bin_cnt_col,
    facet_col_spacing=0.04,
    facet_row_spacing=0.04,
    hover_data=hover_cols,
    hover_name=df_preds.index.name,
    range_x=(domain := (-4, 4) if cli_args.energy_type == Key.each else (-5, 3)),
    range_y=domain,
    category_orders={
        facet_col: sorted([m.label for m in models_to_plot], key=legend_order.index)
    },
    color_continuous_scale="agsunset",
    width=280 * n_cols,
    height=230 * n_rows,
)

# decrease marker size
fig.update_traces(marker=dict(size=2))
# manually set colorbar ticks and labels (needed after log1p transform)
tick_vals = [1, 10, 100, 1000, 10_000]
fig.layout.coloraxis.colorbar.update(
    tickvals=np.log1p(tick_vals), ticktext=list(map(str, tick_vals))
)

# iterate over subplots and set new title
for idx, anno in enumerate(fig.layout.annotations, start=1):
    traces = [t for t in fig.data if t.xaxis == f"x{idx if idx > 1 else ''}"]

    model = anno.text.split("=", 1)[1]
    if model not in df_preds:
        print(f"Unexpected {model=}, not in {[m.label for m in models_to_plot]}")
        continue

    # add MAE and R2 to subplot titles
    MAE = dfs_metrics[test_subset].loc["MAE", model]
    R2 = dfs_metrics[test_subset].loc["R2", model]
    sub_title = f"{model} · {MAE=:.2f} · R<sup>2</sup>={R2:.2f}"
    fig.layout.annotations[idx - 1].text = sub_title

# add transparent rectangle with TN, TP, FN, FP labels in each quadrant
annotate_quadrants = True
if e_true_col == MbdKey.each_true:
    # add dashed quadrant separators
    fig.add_vline(x=0, line=dict(width=0.5, dash="dash"))
    fig.add_hline(y=0, line=dict(width=0.5, dash="dash"))

    if annotate_quadrants:
        for sign_x, sign_y, label, color in (
            (-1, -1, "TP", clf_colors[0]),
            (-1, 1, "FN", clf_colors[1]),
            (1, -1, "FP", clf_colors[2]),
            (1, 1, "TN", clf_colors[3]),
        ):
            # instead of coloring points in each quadrant, we can add a transparent
            # background to each quadrant (looks worse maybe than coloring points)
            fig.add_shape(
                type="rect",
                x0=0,
                y0=0,
                x1=sign_x * 100,
                y1=sign_y * 100,
                fillcolor=color,
                opacity=0.05,
                layer="below",
                row="all",
                col="all",
            )
            fig.add_annotation(
                x=(domain[0] if sign_x < 0 else domain[1]),
                y=(domain[0] if sign_y < 0 else domain[1]),
                xshift=-20 * sign_x,
                yshift=-15 * sign_y if sign_x != sign_y else -70 * sign_y,
                text=label,
                showarrow=False,
                font=dict(size=14, color=color),
                row="all",
                col="all",
            )

pmv.powerups.add_identity_line(fig)

# remove legend title and place legend centered above subplots, increase marker size
fig.layout.legend.update(
    title="", orientation="h", x=0.5, xanchor="center", y=1.15, itemsizing="constant"
)

# Create shared x and y axis titles
axis_titles = dict(xref="paper", yref="paper", showarrow=False, font_size=16)
fig.add_annotation(  # x-axis title
    x=0.5,
    y=0,
    yshift=-50,
    text=e_true_col.label,
    borderpad=5,
    **axis_titles,
)
fig.add_annotation(  # y-axis title
    x=0,
    xshift=-70,
    y=0.5,
    text=e_pred_col.label,
    textangle=-90,
    borderpad=5,
    **axis_titles,
)

# place the colorbar above the subplots
fig.layout.coloraxis.colorbar.update(
    x=0.5, y=1.03, thickness=11, len=0.8, orientation="h"
)

# standardize the margins and template
portrait = n_rows > n_cols
fig.layout.margin.update(l=60, r=10, t=0 if portrait else 10, b=60 if portrait else 10)

axes_kwargs = dict(matches=None, title_text="", showgrid=True, nticks=8)
fig.update_xaxes(**axes_kwargs, range=domain)
fig.update_yaxes(**axes_kwargs, range=domain)
fig.show()


# %%
img_suffix = "" if cli_args.show_non_compliant else "-only-compliant"
img_name = f"{cli_args.energy_type}-parity-models-{n_rows}x{n_cols}{img_suffix}"
pmv.save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf")
