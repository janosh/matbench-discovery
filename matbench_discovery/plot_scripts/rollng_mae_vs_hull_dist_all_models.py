# %%
import pandas as pd

from matbench_discovery import ROOT, today
from matbench_discovery.plot_scripts import df_wbm
from matbench_discovery.plots import plt, rolling_mae_vs_hull_dist

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"


# %%
dfs: dict[str, pd.DataFrame] = {}
dfs["Wren"] = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wren-mp-initial-structures.csv"
).set_index("material_id")
dfs["CGCNN ISRE"] = pd.read_csv(
    # f"{ROOT}/data/2022-06-11-from-rhys/cgcnn-mp-initial-structures.csv"
    f"{ROOT}/models/cgcnn/2022-11-23-test-cgcnn-wbm-IS2RE/cgcnn-ensemble-preds.csv"
).set_index("material_id")
dfs["CGCNN RS2RE"] = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/cgcnn-mp-cse.csv"
).set_index("material_id")
dfs["Voronoi ISRE"] = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/voronoi-mp-initial-structures.csv"
).set_index("material_id")
dfs["Voronoi RS2RE"] = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/voronoi-mp-cse.csv"
).set_index("material_id")
dfs["Wrenformer"] = pd.read_csv(
    f"{ROOT}/models/wrenformer/2022-11-15-wrenformer-IS2RE-preds.csv"
).set_index("material_id")

dfs["megnet"] = (
    pd.read_csv(
        f"{ROOT}/models/megnet/2022-11-18-megnet-wbm-IS2RE/megnet-e-form-preds.csv"
    )
    .set_index("material_id")
    .dropna()
)
dfs["m3gnet"] = pd.read_json(
    f"{ROOT}/models/m3gnet/2022-10-31-m3gnet-wbm-IS2RE.json.gz"
).set_index("material_id")
dfs["bowsr_megnet"] = pd.read_json(
    f"{ROOT}/models/bowsr/2022-11-22-bowsr-megnet-wbm-IS2RE.json.gz"
).set_index("material_id")


# %%
fig, ax = plt.subplots(1, figsize=(10, 9))

target_col = "e_form_per_atom_mp2020_corrected"
e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"

for model_name, df in sorted(dfs.items()):

    if f"e_form_per_atom_{model_name}" in df:
        model_preds = df[f"e_form_per_atom_{model_name}"]
    elif f"{target_col}_pred_ens" in df:
        model_preds = df[f"{target_col}_pred_ens"]
    elif len(pred_cols := df.filter(like=r"_pred_").columns) > 1:
        # make sure we average the expected number of ensemble member predictions
        assert len(pred_cols) == 10, f"{len(pred_cols) = }, expected 10"
        model_preds = df[pred_cols].mean(axis=1)
    elif "e_form_pred" in df:  # voronoi
        model_preds = df.e_form_pred
    else:
        raise ValueError(
            f"No condition matched for {model_name=}, "
            f"which column of {list(df)} holds preds?"
        )
    assert model_preds.isna().sum() < 100

    rolling_mae_vs_hull_dist(
        e_above_hull_pred=model_preds - df_wbm.loc[df.index][target_col],
        e_above_hull_true=df_wbm.loc[df.index][e_above_hull_col],
        ax=ax,
        label=model_name,
    )

# increase line width in legend
legend = ax.legend(frameon=False, loc="lower right")
for line in legend.get_lines():
    line._linewidth *= 3


# %%
img_path = f"{ROOT}/figures/{today}-rolling-mae-vs-hull-dist-compare-models.pdf"
fig.savefig(img_path)
