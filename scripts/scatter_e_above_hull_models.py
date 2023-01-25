# %%
import numpy as np
import plotly.graph_objects as go
from pymatviz.utils import add_identity_line, save_fig
from sklearn.metrics import r2_score

from matbench_discovery import FIGS, STATIC, today
from matbench_discovery.data import PRED_FILENAMES, load_df_wbm_preds
from matbench_discovery.energy import classify_stable
from matbench_discovery.plots import px

__author__ = "Janosh Riebesell"
__date__ = "2022-11-28"


# %%
print(f"loadable models: {list(PRED_FILENAMES)}")
models = sorted(
    "CGCNN, Voronoi Random Forest, Wrenformer, MEGNet, M3GNet, BOWSR MEGNet".split(", ")
)
df_wbm = load_df_wbm_preds(models=models).round(3)

e_form_col = "e_form_per_atom_mp2020_corrected"
e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"
id_col = "material_id"
legend = dict(x=1, y=0, xanchor="right", yanchor="bottom", title=None)


# %%
e_form_preds = "e_form_per_atom_pred"
e_above_hull_preds = "e_above_hull_pred"
var_name = "Model"
hover_cols = (id_col, e_form_col, e_above_hull_col, "formula")

df_melt = df_wbm.melt(
    id_vars=hover_cols,
    value_vars=models,
    var_name=var_name,
    value_name=e_form_preds,
)

df_melt[e_above_hull_preds] = (
    df_melt[e_above_hull_col] - df_melt[e_form_col] + df_melt[e_form_preds]
)


# %%
def _metric_str(xs: list[float], ys: list[float]) -> str:
    # compute MAE and R2 for set of (x, y) pairs
    isna = np.isnan(xs) | np.isnan(ys)
    xs, ys = xs[~isna], ys[~isna]
    MAE = np.abs(xs - ys).mean()
    R2 = r2_score(xs, ys)
    return f"· MAE={MAE:.2f} · R<sup>2</sup>={R2:.2f}"


def _add_metrics_to_legend(fig: go.Figure) -> None:
    for trace in fig.data:
        # initially hide all traces, let users select which models to compare
        trace.visible = "legendonly"
        trace.name = f"{trace.name}{_metric_str(trace.x, trace.y)}"


# %% scatter plot of actual vs predicted e_form_per_atom
fig = px.scatter(
    df_melt.iloc[::5],
    x=e_form_col,
    y=e_form_preds,
    color=var_name,
    hover_data=hover_cols,
    hover_name=id_col,
)

_add_metrics_to_legend(fig)
fig.update_layout(legend=legend)
add_identity_line(fig)
fig.show()


# %%
img_path = f"{FIGS}/{today}-e-form-scatter-models"
# save_fig(fig, f"{img_path}.svelte")


# %% scatter plot of actual vs predicted e_above_hull
fig = px.scatter(
    df_melt.iloc[::5],
    x=e_above_hull_col,
    y=e_above_hull_preds,
    color=var_name,
    hover_data=hover_cols,
    hover_name=id_col,
)

_add_metrics_to_legend(fig)
fig.update_layout(legend=legend)
add_identity_line(fig)
fig.show()


# %%
img_path = f"{FIGS}/{today}-e-above-hull-scatter-models"
# save_fig(fig, f"{img_path}.svelte")


# %% plot all models in separate subplots
true_pos, false_neg, false_pos, true_neg = classify_stable(
    df_melt[e_above_hull_col], df_melt[e_above_hull_preds]
)

df_melt["clf"] = np.array(
    classes := ["true positive", "false negative", "false positive", "true negative"]
)[true_pos * 0 + false_neg * 1 + false_pos * 2 + true_neg * 3]

fig = px.scatter(
    df_melt.iloc[::50],
    x=e_above_hull_col,
    y=e_above_hull_preds,
    facet_col=var_name,
    facet_col_wrap=3,
    facet_col_spacing=0.04,
    facet_row_spacing=0.15,
    hover_data=hover_cols,
    hover_name=id_col,
    color="clf",
    color_discrete_map=dict(zip(classes, ("green", "yellow", "red", "blue"))),
    # opacity=0.4,
    range_x=[-2, 2],
    range_y=[-2, 2],
)

x_title = fig.layout.xaxis.title.text
y_title = fig.layout.yaxis.title.text

# iterate over subplots and set new title
for idx, anno in enumerate(fig.layout.annotations, 1):
    traces = [t for t in fig.data if t.xaxis == f"x{idx if idx > 1 else ''}"]
    xs = np.concatenate([t.x for t in traces])
    ys = np.concatenate([t.y for t in traces])

    model = anno.text.split("=")[1]
    # set new subplot titles (adding MAE and R2)
    fig.layout.annotations[idx - 1].text = f"{model} {_metric_str(xs, ys)}"

    # remove x and y axis titles if not on center row or center column
    fig.layout[f"xaxis{idx}"].title.text = ""
    fig.layout[f"yaxis{idx}"].title.text = ""

    # add vertical and horizontal lines at 0
    fig.add_vline(x=0, line=dict(width=1, dash="dash", color="gray"))
    fig.add_hline(y=0, line=dict(width=1, dash="dash", color="gray"))


fig.update_xaxes(nticks=5)
fig.update_yaxes(nticks=5)

legend = dict(
    title="",  # remove legend title
    itemsizing="constant",  # increase legend marker size
    orientation="h",
    x=0.5,  # place legend centered above subplots
    xanchor="center",
    y=1.2,
    yanchor="top",
)
fig.layout.legend.update(legend)

axis_titles = dict(xref="paper", yref="paper", showarrow=False)
fig.add_annotation(
    x=0.5,
    y=-0.16,
    text=x_title,
    **axis_titles,
)
# add y-axis title
fig.add_annotation(
    x=-0.06,
    y=0.5,
    text=y_title,
    textangle=-90,
    **axis_titles,
)


fig.show()
img_path = f"{STATIC}/{today}-each-scatter-models.webp"
save_fig(fig, img_path, scale=4, width=1000, height=500)
