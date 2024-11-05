"""Analyze symmetry retention across models."""

# %%
import os

import pandas as pd
import plotly.express as px
import pymatviz as pmv
from IPython.display import display
from pymatgen.core import Structure
from pymatviz.enums import Key, Task
from tqdm import tqdm

from matbench_discovery import ROOT, SITE_FIGS, today
from matbench_discovery.data import DataFiles, Model, df_wbm, round_trip_yaml
from matbench_discovery.enums import MbdKey
from matbench_discovery.structure import analyze_symmetry, pred_vs_ref_struct_symmetry

init_spg_col = "init_spg_num"
dft_spg_col = "dft_spg_num"
model_lvl, sym_prop_lvl = "model", "symmetry_property"
df_wbm[init_spg_col] = df_wbm[MbdKey.init_wyckoff].str.split("_").str[2].astype(int)
df_wbm[dft_spg_col] = df_wbm[Key.wyckoff_spglib].str.split("_").str[2].astype(int)
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
    n_structs_for_model = len(dfs_model_structs[model.name])
    print(f"Loaded {n_structs_for_model:,} structures for {model.name}")


# %% Perform symmetry analysis for all model-relaxed structures
dfs_sym_all: dict[str, pd.DataFrame] = {}
df_structs = pd.DataFrame()

for model_name in tqdm(dfs_model_structs, desc="Analyzing model structures"):
    n_structs_for_model = len(dfs_model_structs[model_name].dropna())
    n_structs_analyzed = len(dfs_sym_all.get(model_name, []))
    if n_structs_analyzed / n_structs_for_model > 0.97:
        # skip model if >97% of its structures already analyzed
        continue  # accounts for structures failing symmetry analysis

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
    dfs_sym_all[model_name] = analyze_symmetry(
        df_structs[model_name].dropna().to_dict(),
        pbar=dict(desc=f"Analyzing {model_name} symmetries"),
    )

# Analyze DFT structures
dft_structs = {
    mat_id: Structure.from_dict(cse[Key.structure])
    for mat_id, cse in df_wbm_structs[Key.cse].items()
}
dfs_sym_all[Key.dft] = analyze_symmetry(dft_structs)


# %% Compare symmetry with DFT reference
n_models = len({*dfs_sym_all} - {Key.dft})

for idx, model_name in enumerate({*dfs_sym_all} - {Key.dft}):
    dfs_sym_all[model_name] = pred_vs_ref_struct_symmetry(
        dfs_sym_all[model_name],
        dfs_sym_all[Key.dft],
        df_structs[model_name].dropna().to_dict(),
        dft_structs,
        pbar=dict(desc=f"{idx+1}/{n_models} Comparing DFT vs {model_name} symmetries"),
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
    df_sym_all.xs(Key.n_sym_ops, level=sym_prop_lvl, axis="columns"),
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


# %%
def analyze_symmetry_changes(df_sym_all: pd.DataFrame) -> pd.DataFrame:
    """Analyze how often each model's predicted structure has different symmetry vs DFT.

    Returns:
        pd.DataFrame: DataFrame with columns for fraction of structures where symmetry
            decreased, matched, or increased vs DFT.
    """
    results: dict[str, dict[str, float]] = {}

    for model in df_sym_all.columns.levels[0]:
        if model == Key.dft:  # don't compare DFT to itself
            continue

        spg_diff = df_sym_all[model][MbdKey.spg_num_diff]
        n_sym_ops_diff = df_sym_all[model][MbdKey.n_sym_ops_diff]
        total = len(spg_diff.dropna())

        # Count cases where spacegroup changed
        changed_mask = spg_diff != 0
        # Among changed cases, count whether symmetry increased or decreased
        sym_decreased = (n_sym_ops_diff < 0) & changed_mask
        sym_increased = (n_sym_ops_diff > 0) & changed_mask
        sym_matched = ~changed_mask

        results[model] = {
            str(Key.symmetry_decrease): float(sym_decreased.sum() / total),
            str(Key.symmetry_match): float(sym_matched.sum() / total),
            str(Key.symmetry_increase): float(sym_increased.sum() / total),
        }

    return pd.DataFrame(results).T


def write_symmetry_metrics_to_yaml(model: Model) -> None:
    """Write symmetry metrics to model YAML metadata files."""
    df_rmsd = df_sym_all.xs(
        MbdKey.structure_rmsd_vs_dft, level=sym_prop_lvl, axis="columns"
    )
    if model.name not in df_rmsd:
        print(f"No RMSD column for {model.name}")
        return

    # Calculate RMSD
    rmsd = float(df_rmsd[model.name].mean(axis=0).round(3))

    # Calculate symmetry change statistics
    df_sym_changes = analyze_symmetry_changes(df_sym_all)
    if model.name not in df_sym_changes.index:
        print(f"No symmetry data for {model.name}")
        return

    sym_changes = df_sym_changes.loc[model.name].round(3).to_dict()

    # Combine metrics
    symmetry_metrics = {str(Key.rmsd): rmsd, str(Key.symmetry_change): sym_changes}

    yaml_path = f"{Model.base_dir}/{model.url}"
    with open(yaml_path) as file:  # Load existing metadata
        model_metadata = round_trip_yaml.load(file)

    model_metadata.setdefault("metrics", {})[Task.geo_opt] = symmetry_metrics

    with open(yaml_path, mode="w") as file:  # Write back to file
        round_trip_yaml.dump(model_metadata, file)


# %% Print summary of symmetry changes
df_sym_changes = analyze_symmetry_changes(df_sym_all)
print("\nSymmetry changes vs DFT (as fraction of total structures):")

display(
    df_sym_changes.round(3)
    .rename(columns=lambda col: col.removeprefix("symmetry_"))
    .style.format("{:.1%}")
    .background_gradient(cmap="Oranges")
    .background_gradient(cmap="Greens", subset="match")
)


# %%
if __name__ == "__main__":
    for model in Model:
        write_symmetry_metrics_to_yaml(model)

    # %% plot ML vs DFT relaxed spacegroup correspondence as sankey diagrams
    df_spg = df_sym_all.xs(Key.spg_num, level=sym_prop_lvl, axis="columns")
    for model in {*df_spg} - {Key.dft}:
        # get most common pairs of DFT/Model spacegroups
        common_dft_spgs, common_model_spgs = zip(
            *df_spg[[Key.dft, model]].value_counts().head(10).index
        )
        model_label = getattr(Model, model).label
        df_spg_common = (
            df_spg.query(f"dft in {common_dft_spgs} and {model} in {common_model_spgs}")
            .sort_values(by=Key.dft)
            .rename(columns={model: model_label, Key.dft: Key.dft.label})
        )

        fig = pmv.sankey_from_2_df_cols(
            df_spg_common.reset_index(),
            [Key.dft.label, model_label],
            label_standoff=35,
        )
        fig.show()
        pmv.save_fig(fig, f"{SITE_FIGS}/spg-sankey-{model.replace('_', '-')}.svelte")
