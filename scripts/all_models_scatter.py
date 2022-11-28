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

target_col = "e_form_per_atom_mp2020_corrected"
e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"
id_col = "material_id"


# %%
value_name = "e_form_per_atom_pred"
var_name = "Model"
hover_cols = (id_col, target_col, e_above_hull_col, "formula")
fig = px.scatter(
    df_wbm.query(f"{target_col} < 5")
    .iloc[::5]
    .melt(
        id_vars=hover_cols,
        value_vars=sorted(models),
        var_name=var_name,
        value_name=value_name,
    ),
    x=target_col,
    y=value_name,
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
