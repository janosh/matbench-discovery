# %%
from pymatviz.utils import add_identity_line
from sklearn.metrics import r2_score

from matbench_discovery import FIGS, today
from matbench_discovery.data import PRED_FILENAMES, load_df_wbm_with_preds
from matbench_discovery.plots import px, write_html

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


# %% scatter plot of actual vs predicted e_form_per_atom
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

fig = px.scatter(
    df_melt.iloc[::10],
    x=e_form_col,
    y=e_form_preds,
    color=var_name,
    hover_data=hover_cols,
    hover_name=id_col,
)


for trace in fig.data:
    # initially hide all traces, let users select which models to compare
    trace.visible = "legendonly"
    # add MAE and R2 to legend
    model = trace.name
    MAE = (df_wbm[e_form_col] - df_wbm[model]).abs().mean()
    R2 = r2_score(*df_wbm[[e_form_col, model]].dropna().to_numpy().T)
    trace.name = f"{model} {MAE=:.2} {R2=:.2}"

fig.update_layout(legend=dict(x=1, y=0, xanchor="right", yanchor="bottom"))
add_identity_line(fig)
fig.show()


# %% scatter plot of actual vs predicted e_above_hull
fig = px.scatter(
    df_melt.iloc[::10],
    x=e_above_hull_col,
    y=e_above_hull_preds,
    color=var_name,
    hover_data=hover_cols,
    hover_name=id_col,
)
for trace in fig.data:
    trace.visible = "legendonly"

    model = trace.name
    MAE = (df_wbm[e_above_hull_col] - df_wbm[e_form_col] + df_wbm[model]).abs().mean()
    R2 = r2_score(*df_wbm[[e_above_hull_col, model]].dropna().to_numpy().T)
    trace.name = f"{model} {MAE=:.2} {R2=:.2g}"

fig.update_layout(legend=dict(x=1, y=0, xanchor="right", yanchor="bottom"))
add_identity_line(fig)
fig.show()


# %%
img_path = f"{FIGS}/{today}-rolling-mae-vs-hull-dist-compare-models"
fig.savefig(f"{img_path}.pdf")
write_html(fig, f"{img_path}.svelte")
