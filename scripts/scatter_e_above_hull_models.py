# %%
import numpy as np
import plotly.graph_objects as go
from pymatviz.utils import add_identity_line, save_fig
from sklearn.metrics import r2_score

from matbench_discovery import FIGS, STATIC, today
from matbench_discovery.data import PRED_FILENAMES, load_df_wbm_with_preds
from matbench_discovery.energy import classify_stable
from matbench_discovery.plots import px

__author__ = "Janosh Riebesell"
__date__ = "2022-11-28"


# %%
print(f"loadable models: {list(PRED_FILENAMES)}")
models = sorted(
    "CGCNN, Voronoi RF, Wrenformer, MEGNet, M3GNet, BOWSR MEGNet".split(", ")
)
df_wbm = load_df_wbm_with_preds(models=models).round(3)

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
def _metric_str(model_name: str) -> str:
    model_pred = df_wbm[e_above_hull_col] - (df_wbm[e_form_col] - df_wbm[model_name])
    MAE = (df_wbm[e_above_hull_col] - model_pred).abs().mean()
    isna = df_wbm[e_above_hull_col].isna() | model_pred.isna()
    R2 = r2_score(df_wbm[e_above_hull_col][~isna], model_pred[~isna])
    return f"{model_name} · {MAE=:.2f} · R<sup>2</sup>={R2:.2f}"


def _add_metrics_to_legend(fig: go.Figure) -> None:
    for trace in fig.data:
        # initially hide all traces, let users select which models to compare
        trace.visible = "legendonly"
        # add MAE and R2 to legend
        model = trace.name
        trace.name = _metric_str(model)


# %% scatter plot of actual vs predicted e_form_per_atom
fig = px.scatter(
    df_melt.iloc[::10],
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
# fig.write_image(f"{img_path}.pdf")
save_fig(fig, f"{img_path}.svelte")


# %% scatter plot of actual vs predicted e_above_hull
fig = px.scatter(
    df_melt.iloc[::10],
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
# fig.write_image(f"{img_path}.pdf")
save_fig(fig, f"{img_path}.svelte")


# %% plot all models in separate subplots
true_pos, false_neg, false_pos, true_neg = classify_stable(
    df_melt[e_above_hull_col], df_melt[e_above_hull_preds]
)

df_melt["clf"] = np.array(
    classes := ["true positive", "false negative", "false positive", "true negative"]
)[true_pos * 0 + false_neg * 1 + false_pos * 2 + true_neg * 3]

fig = px.scatter(
    df_melt.iloc[::10],
    x=e_above_hull_col,
    y=e_above_hull_preds,
    facet_col=var_name,
    facet_col_wrap=3,
    hover_data=hover_cols,
    hover_name=id_col,
    color="clf",
    color_discrete_map=dict(zip(classes, ("green", "yellow", "red", "blue"))),
    opacity=0.4,
)

# iterate over subplots and set new title
for idx, model in enumerate(models, 1):
    # find index of annotation belonging to model
    anno_idx = [a.text for a in fig.layout.annotations].index(f"Model={model}")
    assert anno_idx >= 0, f"could not find annotation for {model}"

    # set new subplot titles (adding MAE and R2)
    fig.layout.annotations[anno_idx].text = _metric_str(model)

    # remove x and y axis titles if not on center row or center column
    if idx != 2:
        fig.layout[f"xaxis{idx}"].title.text = ""
    if idx > 1:
        fig.layout[f"yaxis{idx}"].title.text = ""

    # add vertical and horizontal lines at 0
    fig.add_vline(x=0, line=dict(width=1, dash="dash", color="gray"))
    fig.add_hline(y=0, line=dict(width=1, dash="dash", color="gray"))

fig.update_layout(showlegend=False)
fig.update_xaxes(nticks=5)
fig.update_yaxes(nticks=5)

fig.show()
img_path = f"{STATIC}/{today}-each-scatter-models.png"
# save_fig(fig, img_path, scale=4, width=1000, height=500)
