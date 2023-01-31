# %%
import numpy as np
from pymatviz.utils import add_identity_line, save_fig

from matbench_discovery import FIGS, STATIC, today
from matbench_discovery.data import PRED_FILENAMES, load_df_wbm_preds
from matbench_discovery.energy import classify_stable, stable_metrics
from matbench_discovery.plots import clf_color_map, clf_labels, px

__author__ = "Janosh Riebesell"
__date__ = "2022-11-28"


# %%
print(f"loadable models: {list(PRED_FILENAMES)}")
models = sorted(
    "CGCNN, Voronoi Random Forest, Wrenformer, MEGNet, M3GNet, BOWSR MEGNet".split(", ")
)
df_wbm = load_df_wbm_preds(models).round(3)

e_form_col = "e_form_per_atom_mp2020_corrected"
each_true_col = "e_above_hull_mp2020_corrected_ppd_mp"
e_form_pred_col = "e_form_per_atom_pred"
each_pred_col = "e_above_hull_pred"
legend = dict(x=1, y=0, xanchor="right", yanchor="bottom", title=None)


# %%
facet_col = "Model"
hover_cols = (df_wbm.index.name, e_form_col, each_true_col, "formula")

df_melt = df_wbm.melt(
    id_vars=hover_cols,
    value_vars=models,
    var_name=facet_col,
    value_name=e_form_pred_col,
)

df_melt[each_pred_col] = (
    df_melt[each_true_col] + df_melt[e_form_pred_col] - df_melt[e_form_col]
)


# %%
def _metric_str(xs: list[float], ys: list[float]) -> str:
    MAE, R2 = (stable_metrics(xs, ys)[x] for x in ["MAE", "R2"])
    return f"· {MAE=:.2f} · R<sup>2</sup>={R2:.2f}"


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
    trace.name = f"{trace.name}{_metric_str(trace.x, trace.y)}"
fig.update_layout(legend=legend)
add_identity_line(fig)
fig.show()


# %%
img_path = f"{FIGS}/{today}-e-form-scatter-models"
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
    trace.name = f"{trace.name}{_metric_str(trace.x, trace.y)}"
fig.update_layout(legend=legend)
add_identity_line(fig)
fig.show()


# %%
img_path = f"{FIGS}/{today}-e-above-hull-scatter-models"
# save_fig(fig, f"{img_path}.svelte")


# %% plot all models in separate subplots
true_pos, false_neg, false_pos, true_neg = classify_stable(
    df_melt[each_true_col], df_melt[each_pred_col]
)

df_melt[(clf_col := "classified")] = np.array(clf_labels)[
    true_pos * 0 + false_neg * 1 + false_pos * 2 + true_neg * 3
]
xy_max = 1.5

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
    opacity=0.4,
    range_x=(-xy_max, xy_max),
    range_y=(-xy_max, xy_max),
    height=1000,
)


x_title = fig.layout.xaxis.title.text  # used in annotations below
y_title = fig.layout.yaxis.title.text

# iterate over subplots and set new title
for idx, anno in enumerate(fig.layout.annotations, 1):
    traces = [t for t in fig.data if t.xaxis == f"x{idx if idx > 1 else ''}"]
    assert len(traces) == 4, f"Expected 4 traces, got {len(traces)=}"
    xs = np.concatenate([t.x for t in traces])
    ys = np.concatenate([t.y for t in traces])
    metrics = stable_metrics(xs, ys)
    MAE, R2 = metrics["MAE"], metrics["R2"]

    model = anno.text.split("=")[1]
    assert model in models, f"Unexpected {model=} not in {models=}"
    # set new subplot titles (adding MAE and R2)
    fig.layout.annotations[idx - 1].text = f"{model}{_metric_str(xs, ys)}"

    # remove subplot x and y axis titles
    fig.layout[f"xaxis{idx}"].title.text = ""
    fig.layout[f"yaxis{idx}"].title.text = ""

    # add transparent rectangle with TN, TP, FN, FP labels in each quadrant
    for sign_x, sign_y, color, label in (
        (-1, -1, "lightseagreen", "TP"),
        (-1, 1, "lightgoldenrodyellow", "FN"),
        (1, -1, "lightsalmon", "FP"),
        (1, 1, "dodgerblue", "TN"),
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
        #     opacity=0.5,
        #     layer="below",
        #     xref=f"x{idx}",
        #     yref=f"y{idx}",
        # )
        fig.add_annotation(
            xref=f"x{idx}",
            yref=f"y{idx}",
            x=sign_x * xy_max,
            y=sign_y * xy_max,
            xshift=-20 * sign_x,
            yshift=-20 * sign_y,
            text=label,
            showarrow=False,
            font=dict(size=16, color=color),
        )

    # add dashed quadrant separators
    fig.add_vline(x=0, line=dict(width=0.5, dash="dash"))
    fig.add_hline(y=0, line=dict(width=0.5, dash="dash"))

fig.update_xaxes(nticks=5)
fig.update_yaxes(nticks=5)

legend = dict(
    title="",  # remove legend title
    itemsizing="constant",  # increase legend marker size
    orientation="h",
    x=0.5,  # place legend centered above subplots
    xanchor="center",
    y=1.1,
    yanchor="top",
)
fig.layout.legend.update(legend)

# fig.update_layout(yaxis=dict(scaleanchor="x", scaleratio=1))

axis_titles = dict(xref="paper", yref="paper", showarrow=False)
fig.add_annotation(  # x-axis title
    x=0.5,
    y=-0.06,
    text=x_title,
    **axis_titles,
)
fig.add_annotation(  # y-axis title
    x=-0.05,
    y=0.5,
    text=y_title,
    textangle=-90,
    **axis_titles,
)

fig.show()


# %%
img_path = f"{STATIC}/{today}-each-scatter-models.webp"
save_fig(fig, img_path, scale=4, width=800, height=700)
