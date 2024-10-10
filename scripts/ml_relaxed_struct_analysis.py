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
from matbench_discovery.structure import analyze_symmetry, compare_symmetry

init_spg_col = "init_spg_num"
dft_spg_col = "dft_spg_num"
sym_data_col = "symmetry_data"
n_sym_ops_col = "n_sym_ops"
n_rot_ops_col = "n_rot_ops"
n_trans_ops_col = "n_trans_ops"
df_wbm[init_spg_col] = df_wbm[MbdKey.init_wyckoff].str.split("_").str[2].astype(int)
df_wbm[dft_spg_col] = df_wbm["wyckoff_spglib"].str.split("_").str[2].astype(int)
module_dir = os.path.dirname(__file__)

debug_mode: int = 100  # limit the number of structures loaded per model to this number

retained = (df_wbm[init_spg_col] == df_wbm[dft_spg_col]).sum()

print(
    f"DFT init structure symmetry retained: {retained}/{len(df_wbm)} "
    f"({retained/len(df_wbm):.2%})"
)

csv_path = f"{ROOT}/data/2024-10-08-all-models-symmetry-analysis.csv.gz"
if os.path.isfile(csv_path):
    df_sym_all = pd.read_csv(csv_path, header=[0, 1], index_col=0)


# %%
df_wbm_structs = pd.read_json(DataFiles.wbm_computed_structure_entries.path)
df_wbm_structs = df_wbm_structs.set_index(Key.mat_id)
if debug_mode:
    df_wbm_structs = df_wbm_structs.head(debug_mode)


# %% Load all available model-relaxed structures for all models
dfs_model_structs: dict[str, pd.DataFrame] = {}
json_struct_paths = {
    model.name: json_path
    for model in Model
    if os.path.isfile(json_path := model.path.replace(".csv.gz", ".json.gz"))
}

pbar = tqdm(json_struct_paths, desc="Loading model structures")
for model_name in pbar:
    pbar.set_postfix_str(model_name)
    df_model = pd.read_json(json_struct_paths[model_name]).set_index(Key.mat_id)
    if debug_mode:
        df_model = df_model.head(debug_mode)
    dfs_model_structs[model_name] = df_model
    print(f"Loaded {len(dfs_model_structs[model_name]):,} structures for {model_name}")


# %% Perform symmetry analysis for all model-relaxed structures
dfs_sym_all: dict[str, pd.DataFrame] = {}
ml_structs: dict[str, list[Structure]] = {}
for model_name in tqdm(dfs_model_structs, desc="Analyzing model structures"):
    df_model = dfs_model_structs[model_name]
    struct_col = next(col for col in df_model if Key.structure in col)
    ml_structs[model_name] = [
        Structure.from_dict(struct_dict) for struct_dict in df_model[struct_col]
    ]
    dfs_sym_all[model_name] = analyze_symmetry(ml_structs[model_name])

# Analyze DFT structures
dft_structs = [
    Structure.from_dict(cse[Key.structure]) for cse in df_wbm_structs[Key.cse]
]
dfs_sym_all[Key.dft] = analyze_symmetry(dft_structs)


# %% Compare symmetry with DFT reference
for model_name in dfs_sym_all:
    if model_name != Key.dft:
        dfs_sym_all[model_name] = compare_symmetry(
            dfs_sym_all[model_name],
            dfs_sym_all[Key.dft],
            ml_structs[model_name],
            dft_structs,
        )


# %% Combine all dataframes
df_sym_all = pd.concat(
    dfs_sym_all, axis="columns", names=["model", "symmetry_property"]
)
df_sym_all = df_sym_all.convert_dtypes()


# %% Save results to CSV
csv_path = f"{ROOT}/data/{today}-all-models-symmetry-analysis.csv.gz"
df_sym_all.to_csv(csv_path)


# %% Plot violin plot of RMSD vs DFT
fig_rmsd = px.violin(
    df_sym_all.xs(Key.rmsd, level="symmetry_property", axis="columns")
    .round(3)
    .dropna(),
    title="RMSD of ML Force Field Relaxed Structures vs DFT Relaxed Structures",
    orientation="h",
    color="model",
    box=True,
)
fig_rmsd.layout.xaxis.title = "Model"
fig_rmsd.layout.yaxis.title = "RMSD (Å)"
fig_rmsd.layout.margin.t = 50
fig_rmsd.update_traces(orientation="h", side="positive", width=1.8)

# add annotation for mean for each model
for model, srs_rmsd in df_sym_all.xs(
    Key.rmsd, level="symmetry_property", axis="columns"
).items():
    mean_rmsd = srs_rmsd.mean()
    model_color = next(
        trace["marker"]["color"] for trace in fig_rmsd.data if trace["name"] == model
    )
    # percentiles = [0.5, 0.95]
    for label, value, offset in [
        ("Mean", mean_rmsd, 40),
        # ("90th percentile", srs_rmsd.quantile(0.95), 25),
    ]:
        fig_rmsd.add_annotation(
            x=value,
            y=model,
            text=f"{label}: {value:.3f}",
            ax=offset,
            ay=-offset,
            font=dict(size=12, color=model_color),
        )

fig_rmsd.layout.height = 600
fig_rmsd.layout.legend.update(x=1, y=1, xanchor="right", yanchor="top")
fig_rmsd.show()

pmv.save_fig(fig_rmsd, f"{module_dir}/{today}-rmsd-violin.pdf")


# %% calculate number of model spacegroups agreeing with DFT-relaxed spacegroup
dft_spg_diff_col = "dft_spg_diff"
# iterate over model spacegroup columns
for model, spg_col in df_sym_all.xs(
    Key.spg_num, level="symmetry_property", axis="columns"
).items():
    df_sym_all[(model, dft_spg_diff_col)] = spg_col - df_wbm[dft_spg_col]


# %%
avg_spg_diff = (
    df_sym_all.xs(dft_spg_diff_col, level="symmetry_property", axis="columns")
    .mean(axis=0)
    .round(1)
)
print(f"Average spacegroup number difference vs DFT for each model: {avg_spg_diff}")


# %% violin plot of spacegroup number diff vs DFT
fig_sym = px.violin(
    df_sym_all.xs(dft_spg_diff_col, level="symmetry_property", axis="columns"),
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
fig_sym_ops.layout.yaxis.title = "Number of Symmetry Operations"
fig_sym_ops.layout.margin.t = 50
fig_sym_ops.layout.height = 600
fig_sym_ops.update_traces(orientation="h", side="positive", width=1.8)
fig_sym_ops.show()
pmv.save_fig(fig_sym_ops, f"{module_dir}/{today}-sym-ops-violin.pdf")


# %% violin plot of number of symmetry operations in ML-relaxed structures vs DFT
fig_sym_ops_diff = px.violin(
    df_sym_all.drop(Key.dft, level="model", axis="columns").xs(
        "n_sym_ops_diff", level="symmetry_property", axis="columns"
    ),
    title="Difference in number of Symmetry Operations in ML-relaxed Structures vs DFT",
    orientation="h",
    color="model",
)
fig_sym_ops_diff.layout.xaxis.title = "Number of Symmetry Operations"
fig_sym_ops_diff.layout.yaxis.title = "Model"
fig_sym_ops_diff.layout.margin.t = 50
fig_sym_ops_diff.update_traces(orientation="h", side="positive", width=1.8)

fig_sym_ops_diff.show()
pmv.save_fig(fig_sym_ops_diff, f"{module_dir}/{today}-sym-ops-diff-violin.pdf")
