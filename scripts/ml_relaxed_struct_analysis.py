"""Analyze symmetry retention across models."""

# %%
import os

import pandas as pd
import plotly.express as px
import pymatviz as pmv
from pymatgen.core import Structure
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import ROOT, today
from matbench_discovery.data import DataFiles, Model, df_wbm
from matbench_discovery.enums import MbdKey
from matbench_discovery.structure import analyze_symmetry, pred_vs_ref_struct_symmetry

init_spg_col = "init_spg_num"
dft_spg_col = "dft_spg_num"
n_sym_ops_col = "n_sym_ops"
model_lvl, sym_prop_lvl = "model", "symmetry_property"
df_wbm[init_spg_col] = df_wbm[MbdKey.init_wyckoff].str.split("_").str[2].astype(int)
df_wbm[dft_spg_col] = df_wbm["wyckoff_spglib"].str.split("_").str[2].astype(int)
module_dir = os.path.dirname(__file__)

# limit the number of structures loaded per model to this number, 0 for no limit
debug_mode: int = 0

retained = (df_wbm[init_spg_col] == df_wbm[dft_spg_col]).sum()

print(
    f"DFT-relaxed structures that retained initial structure's spacegroup: "
    f"{retained:,} / {len(df_wbm):,} ({retained/len(df_wbm):.2%})"
)

csv_path = f"{ROOT}/data/2024-10-10-all-models-symmetry-analysis.csv.gz"
if os.path.isfile(csv_path):
    df_sym_all = pd.read_csv(csv_path, header=[0, 1], index_col=0)


# %%
df_wbm_structs = pd.read_json(DataFiles.wbm_computed_structure_entries.path)
df_wbm_structs = df_wbm_structs.set_index(Key.mat_id)
if debug_mode:
    df_wbm_structs = df_wbm_structs.head(debug_mode)


# %% Load all available model-relaxed structures for all models
dfs_model_structs: dict[str, pd.DataFrame] = {}

# pbar = tqdm(Model, desc="Loading model structures")
for model in Model:
    json_path = str(model).replace(".csv.gz", ".json.gz")
    if not os.path.isfile(json_path):
        continue

    # pbar.set_postfix_str(model_name)
    if (
        model.name in dfs_model_structs
        # reload df_model if debug_mode changed
        and len(dfs_model_structs[model.name]) == debug_mode
    ):
        continue
    df_model = pd.read_json(json_path).set_index(Key.mat_id)
    if debug_mode:
        df_model = df_model.head(debug_mode)
    dfs_model_structs[model.name] = df_model
    n_structs = len(dfs_model_structs[model.name])
    print(f"Loaded {n_structs:,} structures for {model.name}")


# %% Perform symmetry analysis for all model-relaxed structures
dfs_sym_all: dict[str, pd.DataFrame] = {}
df_structs = pd.DataFrame()

for model_name in tqdm(dfs_model_structs, desc="Analyzing model structures"):
    df_model = dfs_model_structs[model_name]
    try:
        struct_col = next(col for col in df_model if Key.structure in col)
    except StopIteration:
        print(f"No structure column found for {model_name}")
        continue
    df_structs[model_name] = {
        mat_id: Structure.from_dict(struct_dict)
        for mat_id, struct_dict in df_model[struct_col].items()
    }
    dfs_sym_all[model_name] = analyze_symmetry(df_structs[model_name])

# Analyze DFT structures
dft_structs = {
    mat_id: Structure.from_dict(cse[Key.structure])
    for mat_id, cse in df_wbm_structs[Key.cse].items()
}
dfs_sym_all[Key.dft] = analyze_symmetry(dft_structs)


# %% Compare symmetry with DFT reference
for model_name in {*dfs_sym_all} - {Key.dft}:
    dfs_sym_all[model_name] = pred_vs_ref_struct_symmetry(
        dfs_sym_all[model_name],
        dfs_sym_all[Key.dft],
        df_structs[model_name].to_dict(),
        dft_structs,
    )


# %% Combine all dataframes
df_sym_all = pd.concat(dfs_sym_all, axis="columns", names=[model_lvl, sym_prop_lvl])
df_sym_all = df_sym_all.convert_dtypes()


# %% Save results to CSV
csv_path = f"{ROOT}/data/{today}-all-models-symmetry-analysis.csv.gz"
df_sym_all.to_csv(csv_path)


# %% Plot violin plot of RMSD vs DFT
fig_rmsd = px.violin(
    df_sym_all.xs(MbdKey.structure_rmsd_vs_dft, level=sym_prop_lvl, axis="columns")
    .round(3)
    .dropna(),
    orientation="h",
    color="model",
    box=True,
)
title = "RMSD of ML-relaxed structures vs DFT-relaxed structures"
fig_rmsd.layout.title = dict(text=title, x=0.5)
fig_rmsd.layout.xaxis.title = "Model"
fig_rmsd.layout.yaxis.title = "RMSD (Å)"
fig_rmsd.layout.margin.t = 50
fig_rmsd.update_traces(orientation="h", side="positive", width=1.8)

# add annotation for mean for each model
for model, srs_rmsd in df_sym_all.xs(
    MbdKey.structure_rmsd_vs_dft, level=sym_prop_lvl, axis="columns"
).items():
    mean_rmsd = srs_rmsd.mean()
    model_color = next(
        (trace["marker"]["color"] for trace in fig_rmsd.data if trace["name"] == model),
        "black",
    )
    fig_rmsd.add_annotation(
        x=mean_rmsd,
        y=model,
        text=f"Mean: {mean_rmsd:.3f}",
        ax=0,
        ay=-60,
        font=dict(size=12, color=model_color),
        arrowcolor=model_color,
    )

fig_rmsd.layout.height = len(fig_rmsd.data) * 100
fig_rmsd.layout.showlegend = False
fig_rmsd.show()

pmv.save_fig(fig_rmsd, f"{module_dir}/{today}-rmsd-violin.pdf")


# %% calculate number of model spacegroups agreeing with DFT-relaxed spacegroup
avg_spg_diff = (
    df_sym_all.xs(MbdKey.spg_num_diff, level=sym_prop_lvl, axis="columns")
    .mean(axis=0)
    .round(1)
)
print(f"Average spacegroup number difference vs DFT for each model: {avg_spg_diff}")


# %% violin plot of spacegroup number diff vs DFT
fig_sym = px.violin(
    df_sym_all.xs(MbdKey.spg_num_diff, level=sym_prop_lvl, axis="columns"),
    title="Spacegroup Number Diff vs DFT",
    orientation="h",
    color="model",
)
fig_sym.layout.xaxis.title = "Model"
fig_sym.layout.yaxis.title = "Spacegroup Number Diff vs DFT"
fig_sym.layout.margin.t = 50
fig_sym.layout.height = 600
fig_sym.update_traces(orientation="h", side="positive", width=1.8)

fig_sym.show()
pmv.save_fig(fig_sym, f"{module_dir}/{today}-sym-violin.pdf")


# %% violin plot of number of symmetry operations in ML-relaxed structures
fig_sym_ops = px.violin(
    df_sym_all.xs(n_sym_ops_col, level="symmetry_property", axis="columns"),
    title="Number of Symmetry Operations in ML-relaxed Structures",
    orientation="h",
    color="model",
)
fig_sym_ops.layout.xaxis.title = "Model"
fig_sym_ops.layout.yaxis.title = "Number of Symmetry Operations N<sub>sop</sub>"
fig_sym_ops.layout.margin.t = 50
fig_sym_ops.layout.height = 600
fig_sym_ops.update_traces(orientation="h", side="positive", width=1.8)
fig_sym_ops.show()
pmv.save_fig(fig_sym_ops, f"{module_dir}/{today}-sym-ops-violin.pdf")


# %% violin plot of number of symmetry operations in ML-relaxed structures vs DFT
fig_sym_ops_diff = px.violin(
    df_sym_all.drop(Key.dft, level=model_lvl, axis="columns")
    .xs(MbdKey.n_sym_ops_diff, level=sym_prop_lvl, axis="columns")
    .reset_index(),
    orientation="h",
    color="model",
    hover_name=Key.mat_id,
)
title = "Difference in number of symmetry operations in ML vs DFT-relaxed structures"
fig_sym_ops_diff.layout.title = dict(text=title, x=0.5)
fig_sym_ops_diff.layout.xaxis.title = "N<sub>sop,ML</sub> - N<sub>sop,DFT</sub>"
fig_sym_ops_diff.layout.yaxis.title = "Model"
fig_sym_ops_diff.layout.margin.t = 50
fig_sym_ops_diff.update_traces(orientation="h", side="positive", width=1.8)

fig_sym_ops_diff.show()
pmv.save_fig(fig_sym_ops_diff, f"{module_dir}/{today}-sym-ops-diff-violin.pdf")
