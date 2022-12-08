# %%
from pymatviz.utils import add_identity_line

from matbench_discovery import ROOT, today
from matbench_discovery.load_preds import data_paths, load_df_wbm_with_preds
from matbench_discovery.plots import px

__author__ = "Janosh Riebesell"
__date__ = "2022-11-28"


# %%
print(f"loadable models: {list(data_paths)}")
models = (
    "Wren, CGCNN, CGCNN IS2RE, CGCNN RS2RE, Voronoi RF, "
    "Wrenformer, MEGNet, M3GNet, BOWSR MEGNet"
).split(", ")
df_wbm = load_df_wbm_with_preds(models=models).round(3)

e_form_col = "e_form_per_atom_mp2020_corrected"
e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"
id_col = "material_id"


# %% scatter plot of actual vs predicted e_form_per_atom
e_form_preds = "e_form_per_atom_pred"
e_above_hull_preds = "e_above_hull_pred"
var_name = "Model"
hover_cols = (id_col, e_form_col, e_above_hull_col, "formula")

df_melt = df_wbm.query(f"{e_form_col} < 5").melt(
    id_vars=hover_cols,
    value_vars=sorted(models),
    var_name=var_name,
    value_name=e_form_preds,
)

df_melt[e_above_hull_preds] = (
    df_melt[e_above_hull_col] - df_melt[e_form_col] + df_melt[e_form_preds]
)

fig = px.scatter(
    df_melt.iloc[::5],
    x=e_form_col,
    y=e_form_preds,
    color=var_name,
    hover_data=hover_cols,
    hover_name=id_col,
)
for trace in fig.data:
    trace["visible"] = "legendonly"

add_identity_line(fig)
fig.show()


df_wbm.e_form_per_atom_wbm.hist(bins=100)


# %% scatter plot of actual vs predicted e_above_hull
fig = px.scatter(
    df_melt.iloc[::5],
    x=e_above_hull_col,
    y=e_above_hull_preds,
    color=var_name,
    hover_data=hover_cols,
    hover_name=id_col,
)
for trace in fig.data:
    trace["visible"] = "legendonly"

add_identity_line(fig)
fig.show()


# %%
img_path = f"{ROOT}/figures/{today}-rolling-mae-vs-hull-dist-compare-models.pdf"
# fig.savefig(img_path)
