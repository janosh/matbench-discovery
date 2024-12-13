"""Run this script to add/update geometry optimization analysis for new models to a CSV
file containing all models."""

# %%
import gc
import os
from typing import Final

import pandas as pd
from pymatgen.core import Structure
from pymatviz.enums import Key

from matbench_discovery import ROOT, today
from matbench_discovery.data import DataFiles
from matbench_discovery.models import MODEL_METADATA
from matbench_discovery.structure import analyze_symmetry, pred_vs_ref_struct_symmetry

debug_mode: Final[int] = 0
symprec: Final[float] = 1e-5

csv_path = f"{ROOT}/data/2024-11-29-all-models-geo-opt-analysis-{symprec=}.csv.gz"
csv_path = csv_path.replace("e-0", "e-")

if os.path.isfile(csv_path):
    print(f"Loading existing CSV file {csv_path}")
    df_go = pd.read_csv(csv_path, index_col=0, header=[0, 1])
else:
    print(f"Creating new dataframe to be saved as {csv_path}")
    df_go = pd.DataFrame(columns=pd.MultiIndex(levels=[[], []], codes=[[], []]))


# %% Load WBM reference structures
df_wbm_structs = pd.read_json(DataFiles.wbm_computed_structure_entries.path)
df_wbm_structs = df_wbm_structs.set_index(Key.mat_id)
if debug_mode:
    df_wbm_structs = df_wbm_structs.head(debug_mode)


# %% Analyze DFT structures if not already done
dft_structs = locals().get("dft_structs") or {
    mat_id: Structure.from_dict(cse[Key.structure])
    for mat_id, cse in df_wbm_structs[Key.computed_structure_entry].items()
}
if Key.dft.label in df_go:
    dft_analysis = df_go[Key.dft.label]
else:
    dft_analysis = analyze_symmetry(
        dft_structs,
        pbar=dict(desc=f"Getting DFT symmetries {symprec=}"),
        symprec=symprec,
    )
    for col in dft_analysis:
        df_go[(Key.dft.label, col)] = dft_analysis[col]


# %% Process each model sequentially
for idx, (model_label, model_metadata) in enumerate(MODEL_METADATA.items()):
    prog_str = f"{idx + 1}/{len(MODEL_METADATA)}:"
    if model_label in df_go:
        print(f"{prog_str} {model_label} already analyzed")
        continue

    geo_opt_metrics = model_metadata.get("metrics", {}).get("geo_opt", {})

    # skip models that don't support geometry optimization
    if geo_opt_metrics in ("not applicable", "not available"):
        continue

    ml_relaxed_structs_path = f"{ROOT}/{geo_opt_metrics.get('pred_file')}"
    if not ml_relaxed_structs_path or not os.path.isfile(ml_relaxed_structs_path):
        print(f"⚠️ {model_label}-relaxed structures not found")
        continue

    # Load model structures
    df_model = pd.read_json(ml_relaxed_structs_path).set_index(Key.mat_id)
    if debug_mode:
        df_model = df_model.head(debug_mode)

    try:
        struct_col = next(col for col in df_model if Key.structure in col)
    except StopIteration:
        print(f"⚠️ No structure column found for {model_label}")
        continue

    # Convert structures
    model_structs = {
        mat_id: Structure.from_dict(struct_dict)
        for mat_id, struct_dict in df_model[struct_col].items()
    }

    # Analyze symmetry
    model_analysis = analyze_symmetry(
        model_structs,
        pbar=dict(desc=f"{prog_str} Analyzing {model_label} symmetries"),
        symprec=symprec,
    )

    # Add model results with proper multi-index columns
    for col in model_analysis:
        df_go[(model_label, col)] = model_analysis[col]

    # Compare with DFT reference
    df_sym_change = pred_vs_ref_struct_symmetry(
        model_analysis,
        dft_analysis,
        model_structs,
        dft_structs,
        pbar=dict(desc=f"{prog_str} Comparing DFT vs {model_label} symmetries"),
    )

    # Add comparison results with proper multi-index columns
    for col in df_sym_change:
        df_go[(model_label, col)] = df_sym_change[col]

    # Save after each model is processed
    csv_path = f"{ROOT}/data/{today}-all-models-geo-opt-analysis-{symprec=}.csv.gz"
    df_go.convert_dtypes().to_csv(csv_path)
    print(f"{prog_str} Completed {model_label} and saved results")

    # Free up memory (maybe helps with Jupyter kernel crashes)
    del model_analysis, df_sym_change
    gc.collect()

print("All models processed!")
