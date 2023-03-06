"""Scatter plot of actual vs predicted e_above_hull and e_form_per_atom for all
models. First 2 plots put all models in single figure with togglable traces.
Last plot is split into 2x3 subplots, one for each model.
"""


# %%
import numpy as np
from pymatviz.utils import add_identity_line, save_fig

from matbench_discovery import FIGS, ROOT, STATIC
from matbench_discovery.metrics import classify_stable
from matbench_discovery.plots import clf_color_map, clf_colors, clf_labels, px
from matbench_discovery.preds import (
    df_metrics,
    df_wbm,
    e_form_col,
    each_pred_col,
    each_true_col,
)

__author__ = "Janosh Riebesell"
__date__ = "2022-11-28"

e_form_pred_col = "e_form_per_atom_pred"
legend = dict(x=1, y=0, xanchor="right", yanchor="bottom", title=None)


# %%
facet_col = "Model"
hover_cols = (df_wbm.index.name, e_form_col, each_true_col, "formula")

df_melt = df_wbm.melt(
    id_vars=hover_cols,
    var_name=facet_col,
    value_vars=df_metrics.T.MAE.nsmallest(6).index,  # top 6 models by MAE
    # value_vars=list(df_metrics),  # all models
    value_name=e_form_pred_col,
)

df_melt[each_pred_col] = (
    df_melt[each_true_col] + df_melt[e_form_pred_col] - df_melt[e_form_col]
)


# %% scatter plot of actual vs predicted e_form_per_atom
fig = px.scatter(
    df_melt.iloc[::5],
    x=e_form_col,
    y=e_form_pred_col,
    color=facet_col,
    hover_data=hover_cols,
    hover_name=df_wbm.index.name,
)

for trace in fig.data:
    # initially hide all traces, let users select which models to compare
    trace.visible = "legendonly"
    model = trace.name
    assert model in df_wbm, f"Unexpected {model=} not in {list(df_wbm)=}"
    MAE, R2 = df_metrics[model][["MAE", "R2"]]
    trace.text = f"{model} · {MAE=:.2f} · R<sup>2</sup>={R2:.2f}"

fig.update_layout(legend=legend)
add_identity_line(fig)
fig.show()


# %%
img_path = f"{FIGS}/e-form-scatter-models"
# save_fig(fig, f"{img_path}.svelte")


# %% scatter plot of actual vs predicted e_above_hull
fig = px.scatter(
    df_melt.iloc[::5],
    x=each_true_col,
    y=each_pred_col,
    color=facet_col,
    hover_data=hover_cols,
    hover_name=df_wbm.index.name,
)

for trace in fig.data:
    trace.visible = "legendonly"
    model = trace.name
    assert model in df_wbm, f"Unexpected {model=} not in {list(df_wbm)=}"
    MAE, R2 = df_metrics[model][["MAE", "R2"]]
    trace.text = f"{model} · {MAE=:.2f} · R<sup>2</sup>={R2:.2f}"

fig.update_layout(legend=legend)
add_identity_line(fig)
fig.show()


# %%
img_path = f"{FIGS}/e-above-hull-scatter-models"
# save_fig(fig, f"{img_path}.svelte")


# %% plot all models in separate subplots
true_pos, false_neg, false_pos, true_neg = classify_stable(
    df_melt[each_true_col], df_melt[each_pred_col]
)

df_melt[(clf_col := "classified")] = np.array(clf_labels)[
    true_pos * 0 + false_neg * 1 + false_pos * 2 + true_neg * 3
]
xy_max = 1.6

fig = px.scatter(
    df_melt.iloc[::50],
    x=each_true_col,
    y=each_pred_col,
    facet_col=facet_col,
    facet_col_wrap=2,
    facet_col_spacing=0.02,
    facet_row_spacing=0.04,
    hover_data=hover_cols,
    hover_name=df_wbm.index.name,
    color=clf_col,
    color_discrete_map=clf_color_map,
    # opacity=0.4,
    range_x=(-xy_max, xy_max),
    range_y=(-xy_max, xy_max),
)


x_title = fig.layout.xaxis.title.text  # used in annotations below
y_title = fig.layout.yaxis.title.text

# iterate over subplots and set new title
for idx, anno in enumerate(fig.layout.annotations, 1):
    traces = [t for t in fig.data if t.xaxis == f"x{idx if idx > 1 else ''}"]
    assert len(traces) in (0, 4), f"Plots be empty or have 4 traces, got {len(traces)=}"

    model = anno.text.split("=", 1)[1]
    assert model in df_wbm, f"Unexpected {model=} not in {list(df_wbm)=}"
    # add MAE and R2 to subplot titles
    MAE, R2 = df_metrics[model][["MAE", "R2"]]
    fig.layout.annotations[
        idx - 1
    ].text = f"{model} · {MAE=:.2f} · R<sup>2</sup>={R2:.2f}"

    # remove subplot x and y axis titles
    fig.layout[f"xaxis{idx}"].title.text = ""
    fig.layout[f"yaxis{idx}"].title.text = ""

# add transparent rectangle with TN, TP, FN, FP labels in each quadrant
for sign_x, sign_y, color, label in zip(
    [-1, -1, 1, 1], [-1, 1, -1, 1], clf_colors, ("TP", "FN", "FP", "TN")
):
    # instead of coloring points in each quadrant, we can add a transparent
    # background to each quadrant (looks worse maybe than coloring points)
    # fig.add_shape(
    #     type="rect",
    #     x0=0,
    #     y0=0,
    #     x1=sign_x * 100,
    #     y1=sign_y * 100,
    #     fillcolor=color,
    #     opacity=0.2,
    #     layer="below",
    #     row="all",
    #     col="all",
    # )
    fig.add_annotation(
        x=sign_x * xy_max,
        y=sign_y * xy_max,
        xshift=-20 * sign_x,
        yshift=-20 * sign_y,
        text=label,
        showarrow=False,
        font=dict(size=16, color=color),
        row="all",
        col="all",
    )

# add dashed quadrant separators
fig.add_vline(x=0, line=dict(width=0.5, dash="dash"))
fig.add_hline(y=0, line=dict(width=0.5, dash="dash"))

fig.update_xaxes(nticks=5)
fig.update_yaxes(nticks=5)

# remove legend title and place legend centered above subplots, increase marker size
fig.layout.legend.update(
    title="", orientation="h", x=0.5, xanchor="center", y=1.1, itemsizing="constant"
)

# fig.update_layout(yaxis=dict(scaleanchor="x", scaleratio=1))

axis_titles = dict(xref="paper", yref="paper", showarrow=False)
fig.add_annotation(  # x-axis title
    x=0.5,
    y=-0.1,
    text=x_title,
    **axis_titles,
)
fig.add_annotation(  # y-axis title
    x=-0.07,
    y=0.5,
    text=y_title,
    textangle=-90,
    **axis_titles,
)
fig.update_layout(margin=dict(l=30, r=0, t=0, b=70))

fig.show()


# %%
img_name = "each-scatter-models"
save_fig(fig, f"{STATIC}/{img_name}.webp", scale=4, width=700, height=800)
save_fig(fig, f"{ROOT}/tmp/figures/{img_name}.pdf", width=700, height=800)
