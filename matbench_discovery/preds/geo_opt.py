"""Functions to calculate and save geometry optimization metrics."""

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
csv_path = f"{ROOT}/data/2024-11-07-all-models-symmetry-analysis.csv.gz"


# %%
if os.path.isfile(csv_path):
    df_sym = pd.read_csv(csv_path, header=[0, 1], index_col=0)
else:
    df_sym = pd.DataFrame()
    df_wbm_structs = pd.read_json(DataFiles.wbm_computed_structure_entries.path)
    df_wbm_structs = df_wbm_structs.set_index(Key.mat_id)
    if debug_mode:
        df_wbm_structs = df_wbm_structs.head(debug_mode)

    # %% Load all available model-relaxed structures for all models
    dfs_model_structs: dict[str, pd.DataFrame] = {}

    for model_name, model_metadata in MODEL_METADATA.items():
        if model_name in df_sym:
            print(f"- {model_name} already analyzed")
            continue
        geo_opt_metrics = model_metadata.get("metrics", {}).get("geo_opt", {})
        if geo_opt_metrics in ("not applicable", "not available"):
            continue
        ml_relaxed_structs_path = f"{ROOT}/{geo_opt_metrics.get('pred_file')}"
        if not ml_relaxed_structs_path:
            continue
        if not os.path.isfile(ml_relaxed_structs_path):
            print(f"⚠️ {model_name}-relaxed structures not found")
            continue

        if (
            model_name in dfs_model_structs
            # reload df_model if debug_mode changed
            and (len(dfs_model_structs[model_name]) == debug_mode or debug_mode == 0)
        ):
            continue
        df_model = pd.read_json(ml_relaxed_structs_path).set_index(Key.mat_id)
        if debug_mode:
            df_model = df_model.head(debug_mode)
        dfs_model_structs[model_name] = df_model
        n_structs_for_model = len(dfs_model_structs[model_name])
        print(f"+ Loaded {n_structs_for_model:,} structures for {model_name}")

    # %% Perform symmetry analysis for all model-relaxed structures
    dfs_sym_all: dict[str, pd.DataFrame] = {}
    df_structs = pd.DataFrame()

    for model_name, df_model in dfs_model_structs.items():
        n_structs_for_model = len(df_model.dropna())
        n_structs_analyzed = len(dfs_sym_all.get(model_name, []))
        if n_structs_analyzed / n_structs_for_model > 0.97:
            # skip model if >97% of its structures already analyzed
            continue  # accounts for structures failing symmetry analysis

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

    # %% Analyze DFT structures
    dft_structs = {
        mat_id: Structure.from_dict(cse[Key.structure])
        for mat_id, cse in df_wbm_structs[Key.cse].items()
    }
    if Key.dft.label not in df_sym:
        dfs_sym_all[Key.dft.label] = analyze_symmetry(dft_structs)

    # %% Compare symmetry with DFT reference
    n_models = len({*dfs_sym_all} - {Key.dft})

    for idx, model_name in enumerate({*dfs_sym_all} - {Key.dft}):
        dfs_sym_all[model_name] = pred_vs_ref_struct_symmetry(
            dfs_sym_all[model_name],
            df_sym[Key.dft.label],
            df_structs[model_name].dropna().to_dict(),
            dft_structs,
            pbar=dict(
                desc=f"{idx+1}/{n_models} Comparing DFT vs {model_name} symmetries"
            ),
        )

    # %% Combine all dataframes
    df_sym = df_sym.join(
        pd.concat(dfs_sym_all, axis="columns", names=[Key.model, MbdKey.sym_prop])
    )
    df_sym = df_sym.convert_dtypes()

    # %% Save results to CSV
    csv_path = f"{ROOT}/data/{today}-all-models-symmetry-analysis.csv.gz"
    df_sym.to_csv(csv_path)
