"""Run this script to add/update geometry optimization analysis for new models to a CSV
file containing all models."""

# %%
import os

import pandas as pd
from pymatgen.core import Structure
from pymatviz.enums import Key

from matbench_discovery import ROOT, today
from matbench_discovery.data import DataFiles
from matbench_discovery.enums import MbdKey
from matbench_discovery.models import MODEL_METADATA
from matbench_discovery.structure import analyze_symmetry, pred_vs_ref_struct_symmetry

debug_mode: int = 0
symprec: float = 1e-2
csv_path = f"{ROOT}/data/2024-11-26-all-models-geo-opt-analysis-{symprec=}.csv.gz"
if os.path.isfile(csv_path):
    df_go = pd.read_csv(csv_path, header=[0, 1], index_col=0)
else:
    df_go = pd.DataFrame()


# %%
df_wbm_structs = pd.read_json(DataFiles.wbm_computed_structure_entries.path)
df_wbm_structs = df_wbm_structs.set_index(Key.mat_id)
if debug_mode:
    df_wbm_structs = df_wbm_structs.head(debug_mode)


# %% Load all available model-relaxed structures for all models
dfs_model_structs: dict[str, pd.DataFrame] = {}

for model_label, model_metadata in MODEL_METADATA.items():
    if model_label in df_go:
        print(f"- {model_label} already analyzed")
        continue
    geo_opt_metrics = model_metadata.get("metrics", {}).get("geo_opt", {})

    # skip models that don't support geometry optimization (energy-only models)
    if geo_opt_metrics in ("not applicable", "not available"):
        continue

    ml_relaxed_structs_path = f"{ROOT}/{geo_opt_metrics.get('pred_file')}"
    if not ml_relaxed_structs_path:
        continue
    if not os.path.isfile(ml_relaxed_structs_path):
        print(f"⚠️ {model_label}-relaxed structures not found")
        continue

    if (
        model_label in dfs_model_structs
        # reload df_model if debug_mode changed
        and (len(dfs_model_structs[model_label]) == debug_mode or debug_mode == 0)
    ):
        continue

    df_model = pd.read_json(ml_relaxed_structs_path).set_index(Key.mat_id)
    if debug_mode:
        df_model = df_model.head(debug_mode)
    dfs_model_structs[model_label] = df_model
    n_structs_for_model = len(dfs_model_structs[model_label])
    print(
        f"+ Loaded {n_structs_for_model:,} {model_label} structures from\n"
        f"{ml_relaxed_structs_path}"
    )


# %% Perform symmetry analysis for all model-relaxed structures
dfs_sym_all: dict[str, pd.DataFrame] = {}
df_structs = pd.DataFrame()

for model_label, df_model in dfs_model_structs.items():
    n_structs_for_model = len(df_model.dropna())
    n_structs_analyzed = len(dfs_sym_all.get(model_label, []))
    if n_structs_analyzed / n_structs_for_model > 0.97:
        # skip model if >97% of its structures already analyzed
        continue  # accounts for structures failing symmetry analysis

    try:
        struct_col = next(col for col in df_model if Key.structure in col)
    except StopIteration:
        print(f"No structure column found for {model_label}")
        continue
    df_structs[model_label] = {
        mat_id: Structure.from_dict(struct_dict)
        for mat_id, struct_dict in df_model[struct_col].items()
    }
    dfs_sym_all[model_label] = analyze_symmetry(
        df_structs[model_label].dropna().to_dict(),
        pbar=dict(desc=f"Analyzing {model_label} symmetries"),
        symprec=symprec,
    )


# %% Analyze DFT structures
dft_structs = {
    mat_id: Structure.from_dict(cse[Key.structure])
    for mat_id, cse in df_wbm_structs[Key.cse].items()
}
if Key.dft.label not in df_go:
    dfs_sym_all[Key.dft.label] = analyze_symmetry(dft_structs)


# %% Compare symmetry with DFT reference
n_models = len({*dfs_sym_all} - {Key.dft})

for idx, model_label in enumerate({*dfs_sym_all} - {Key.dft}):
    dfs_sym_all[model_label] = pred_vs_ref_struct_symmetry(
        dfs_sym_all[model_label],
        df_go[Key.dft.label],
        df_structs[model_label].dropna().to_dict(),
        dft_structs,
        pbar=dict(desc=f"{idx+1}/{n_models} Comparing DFT vs {model_label} symmetries"),
    )


# %% Combine existing dataframe with new results and save to CSV
csv_path = f"{ROOT}/data/{today}-all-models-geo-opt-analysis-{symprec=}.csv.gz"
df_go.join(
    pd.concat(dfs_sym_all, axis="columns", names=[Key.model, MbdKey.sym_prop])
).convert_dtypes().to_csv(csv_path)
