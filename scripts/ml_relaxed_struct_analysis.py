"""Analyze symmetry retention across models."""

# %%
import os
import warnings

import numpy as np
import pandas as pd
import plotly.express as px
import pymatviz as pmv
import spglib
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import ROOT, today
from matbench_discovery.data import DataFiles, Model, df_wbm
from matbench_discovery.enums import MbdKey

init_spg_col = "init_spg_num"
dft_spg_col = "dft_spg_num"
sym_data_col = "symmetry_data"
n_sym_ops_col = "n_sym_ops"
n_rot_ops_col = "n_rot_ops"
n_trans_ops_col = "n_trans_ops"
df_wbm[init_spg_col] = df_wbm[MbdKey.init_wyckoff].str.split("_").str[2].astype(int)
df_wbm[dft_spg_col] = df_wbm["wyckoff_spglib"].str.split("_").str[2].astype(int)
module_dir = os.path.dirname(__file__)

debug_mode: int = 0  # limit the number of structures loaded per model to this number

retained = (df_wbm[init_spg_col] == df_wbm[dft_spg_col]).sum()

print(
    f"DFT init structure symmetry retained: {retained}/{len(df_wbm)} "
    f"({retained/len(df_wbm):.2%})"
)

csv_path = f"{ROOT}/data/2024-10-08-all-models-symmetry-analysis.csv.gz"
df_sym_all = pd.read_csv(csv_path, header=[0, 1], index_col=0)

for model in Model:
    struct_json_path = f"{model.path.removesuffix(".csv.gz")}.json.gz"
    print(f"{os.path.isfile(struct_json_path)=}")


# %%
df_wbm_structs = pd.read_json(DataFiles.wbm_computed_structure_entries.path)
df_wbm_structs = df_wbm_structs.set_index(Key.mat_id)
if debug_mode:
    df_wbm_structs = df_wbm_structs.head(debug_mode)


# %% Load all available model-relaxed structures for all models
dfs_model_structs: dict[str, pd.DataFrame] = locals().get("dfs_model_structs", {})
json_struct_paths = {
    model.name: json_path
    for model in Model
    if os.path.isfile(json_path := model.path.replace(".csv.gz", ".json.gz"))
}

pbar = tqdm(json_struct_paths, desc="Loading model structures")
for model_name in pbar:
    pbar.set_postfix_str(model_name)
    if model_name in df_sym_all:
        continue
    df_model = pd.read_json(json_struct_paths[model_name]).set_index(Key.mat_id)
    if debug_mode:
        df_model = df_model.head(debug_mode)
    dfs_model_structs[model_name] = df_model
    print(f"Loaded {len(dfs_model_structs[model_name]):,} structures for {model_name}")


# %% Perform spacegroup analysis and RMSD calculation for all model-relaxed structures
sym_data_col = "symmetry_data"
rmsd_col, max_pair_dist_col = "rmsd_vs_dft", "max_pair_dist"
structure_matcher = StructureMatcher()

dfs_all_structs = {**dfs_model_structs, Key.dft: df_wbm_structs}

# don't iterate over dfs_all_structs itself here because we need to pop keys from it
# which would raise RuntimeError: dictionary changed size during iteration
for model_name in list(dfs_all_structs):
    df_model = dfs_all_structs[model_name]
    try:  # find structure col
        structure_col = next(col for col in df_model if "structure" in col)
    except StopIteration:
        print(f"No structure column found for {model_name}, skipping")
        dfs_all_structs.pop(model_name)
        continue

    if sym_data_col not in df_model:
        df_model[sym_data_col] = None
    if rmsd_col not in df_model:
        df_model[rmsd_col] = np.nan

    for mat_id, struct_dict in tqdm(
        df_model[structure_col].items(),
        desc=f"Processing {model_name} structures",
        total=len(df_model),
    ):
        if not pd.isna(df_model.loc[mat_id, sym_data_col]) and not pd.isna(
            df_model.loc[mat_id, rmsd_col]
        ):
            continue

        model_struct = Structure.from_dict(
            struct_dict
            if model_name != Key.dft
            else df_wbm_structs.loc[mat_id, Key.cse][Key.structure]
        )

        # Symmetry analysis
        cell = (
            model_struct.lattice.matrix,
            model_struct.frac_coords,
            model_struct.atomic_numbers,
        )
        with warnings.catch_warnings():
            warnings.simplefilter(action="ignore", category=spglib.spglib.SpglibError)
            sym_data: spglib.SpglibDataset = spglib.get_symmetry_dataset(cell)
        df_model.at[mat_id, sym_data_col] = sym_data  # noqa: PD008
        if model_name != Key.dft and (dft_cse := df_wbm_structs[Key.cse].get(mat_id)):
            dft_struct = Structure.from_dict(dft_cse[Key.structure])
            rmsd = structure_matcher.get_rms_dist(model_struct, dft_struct)
            df_model.loc[mat_id, (rmsd_col, max_pair_dist_col)] = rmsd


# %%
dfs_sym_all: dict[str, pd.DataFrame] = locals().get("dfs_sym", {})
sym_key_map = {
    "number": Key.spg_num,
    "hall_number": "hall_num",
    "international": "international_spg_name",
    "hall": "hall_symbol",
    "choice": "choice_symbol",
    "pointgroup": Key.point_group,
    "wyckoffs": "wyckoff_symbols",
}
for idx, (model_name, df_model) in enumerate(dfs_all_structs.items(), start=1):
    # TODO figure out how to make this skip criterion depend on df_sym_all instead of
    # dfs_sym
    if dfs_sym_all[model_name].size > 0:
        continue  # skip already processed
    print(f"{idx}/{len(dfs_model_structs)}: processing {model_name} symmetries")
    dfs_sym_all[model_name] = pd.json_normalize(
        df_model[sym_data_col].map(
            lambda sym_data: {
                new_key: getattr(sym_data, old_key)
                for old_key, new_key in sym_key_map.items()
            }
            if sym_data
            else None
        )
    ).convert_dtypes()

    dfs_sym_all[model_name].index = df_model.index

    for mat_id, sym_data in df_model[sym_data_col].items():
        n_rot_ops = len(sym_data.rotations)
        n_trans_ops = len(sym_data.translations)
        dfs_sym_all[model_name].loc[mat_id, n_sym_ops_col] = n_rot_ops + n_trans_ops
        dfs_sym_all[model_name].loc[mat_id, n_rot_ops_col] = n_rot_ops
        dfs_sym_all[model_name].loc[mat_id, n_trans_ops_col] = n_trans_ops


model_lvl, sym_prop_lvl = "model", "symmetry_property"
df_sym_all = pd.concat(
    dfs_sym_all, axis="columns", names=[model_lvl, sym_prop_lvl]
).convert_dtypes()


# %% Add RMSD data to the main DataFrame
for model_name, df_model in dfs_model_structs.items():
    df_sym_all[(model_name, rmsd_col)] = df_model[rmsd_col]

df_sym_all = df_sym_all.convert_dtypes()


# %% subtract DFT number of symmetry operations from ML number of symmetry operations
n_sym_ops_diff_col = "n_sym_ops_diff"
for model, n_sym_ops_srs in df_sym_all.xs(
    n_sym_ops_col, level=sym_prop_lvl, axis="columns"
).items():
    df_sym_all[(model, n_sym_ops_diff_col)] = (
        n_sym_ops_srs - dfs_sym_all[Key.dft][n_sym_ops_col]
    )


# %% Save results to CSV
csv_path = f"{ROOT}/data/{today}-all-models-symmetry-analysis.csv.gz"
df_sym_all.to_csv(csv_path)

# df_sym_all = pd.read_csv(csv_path, header=[0, 1], index_col=0)


# %% Plot violin plot of RMSD vs DFT
fig_rmsd = px.violin(
    df_sym_all.xs(rmsd_col, level=sym_prop_lvl, axis="columns").round(3).dropna(),
    title="RMSD of ML Force Field Relaxed Structures vs DFT Relaxed Structures",
    orientation="h",
    color=model_lvl,
    box=True,
)
fig_rmsd.layout.xaxis.title = "Model"
fig_rmsd.layout.yaxis.title = "RMSD (Å)"
fig_rmsd.layout.margin.t = 50
fig_rmsd.update_traces(orientation="h", side="positive", width=1.8)

# add annotation for mean for each model
for model, srs in df_sym_all.xs(rmsd_col, level=sym_prop_lvl, axis="columns").items():
    mean_rmsd = srs.mean()
    model_color = next(
        trace["marker"]["color"] for trace in fig_rmsd.data if trace["name"] == model
    )
    # percentiles = [0.5, 0.95]
    for label, value, offset in [
        ("Mean", mean_rmsd, 40),
        # ("90th percentile", srs.quantile(0.95), 25),
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
    Key.spg_num, level=sym_prop_lvl, axis="columns"
).items():
    df_sym_all[(model, dft_spg_diff_col)] = spg_col - df_wbm[dft_spg_col]


# %%
avg_spg_diff = (
    df_sym_all.xs(dft_spg_diff_col, level=sym_prop_lvl, axis="columns")
    .mean(axis=0)
    .round(1)
)
print(f"Average spacegroup number difference vs DFT for each model: {avg_spg_diff}")


# %% violin plot of spacegroup number diff vs DFT
fig_sym = px.violin(
    df_sym_all.xs(dft_spg_diff_col, level=sym_prop_lvl, axis="columns"),
    title="Spacegroup Number Diff vs DFT",
    orientation="h",
    color=model_lvl,
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
    df_sym_all.xs(n_sym_ops_col, level=sym_prop_lvl, axis="columns"),
    title="Number of Symmetry Operations in ML-relaxed Structures",
    orientation="h",
    color=model_lvl,
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
    df_sym_all.drop(Key.dft, level=model_lvl, axis="columns").xs(
        n_sym_ops_diff_col, level=sym_prop_lvl, axis="columns"
    ),
    title="Difference in number of Symmetry Operations in ML-relaxed Structures vs DFT",
    orientation="h",
    color=model_lvl,
)
fig_sym_ops_diff.layout.xaxis.title = "Number of Symmetry Operations"
fig_sym_ops_diff.layout.yaxis.title = "Model"
fig_sym_ops_diff.layout.margin.t = 50
fig_sym_ops_diff.update_traces(orientation="h", side="positive", width=1.8)

fig_sym_ops_diff.show()
pmv.save_fig(fig_sym_ops_diff, f"{module_dir}/{today}-sym-ops-diff-violin.pdf")
