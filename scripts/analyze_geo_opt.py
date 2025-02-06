"""Run this script to add/update geometry optimization analysis for new models to
individual CSV files (one per symprec value) in a model's directory.

Output files will have the same name as the input file containing the model's relaxed
structures, but with the symprec value appended to the filename.
"""

# %%
import importlib
import importlib.metadata
import os
from typing import Any, Final

import pandas as pd
from pymatgen.core import Structure
from pymatviz.enums import Key

from matbench_discovery import ROOT
from matbench_discovery.data import DataFiles, Model
from matbench_discovery.metrics import geo_opt
from matbench_discovery.models import MODEL_METADATA
from matbench_discovery.structure import analyze_symmetry, pred_vs_ref_struct_symmetry

# activate debug mode by setting to any number > 0, only that many structures will be
# analyzed
debug_mode: Final[int] = 0
symprec: Final[float] = 1e-2
symprec_str = f"symprec={symprec:.0e}".replace("e-0", "e-")

moyo_version = f"moyo={importlib.metadata.version('moyopy')}"
dft_csv_path = f"{ROOT}/data/wbm/dft-geo-opt-{symprec_str}-{moyo_version}.csv.gz"


# %% Load WBM reference structures (this takes a while)
wbm_cse_path = DataFiles.wbm_computed_structure_entries.path
df_wbm_structs: pd.DataFrame = pd.read_json(wbm_cse_path).set_index(Key.mat_id)
dft_structs: dict[str, Structure] = {
    mat_id: Structure.from_dict(cse[Key.structure])
    for mat_id, cse in df_wbm_structs[Key.computed_structure_entry].items()
}

if debug_mode:
    df_wbm_structs = df_wbm_structs.head(debug_mode)


# %%
if os.path.isfile(dft_csv_path):
    dft_analysis = pd.read_csv(dft_csv_path)
else:
    dft_analysis = analyze_symmetry(
        dft_structs,
        pbar=dict(desc=f"Getting DFT symmetries {symprec=}"),
        symprec=symprec,
    )
    dft_analysis.to_csv(dft_csv_path)


# %% Process each model sequentially
for idx, (model_label, model_metadata) in enumerate(MODEL_METADATA.items()):
    prog_str = f"{idx + 1}/{len(MODEL_METADATA)}:"
    model = Model.from_label(model_label)

    geo_opt_metrics: dict[str, Any] = model_metadata.get("metrics", {}).get(
        "geo_opt", {}
    )

    # skip models that don't support geometry optimization
    if geo_opt_metrics in ("not applicable", "not available"):
        continue

    if not model.geo_opt_path:
        print(f"⚠️ {model_label} has no relaxed structures file")
        continue

    if not os.path.isfile(ml_relaxed_structs_path := model.geo_opt_path):
        print(
            f"⚠️ {model_label}-relaxed structures not found, expected "
            f"at {ml_relaxed_structs_path}"
        )
        continue

    # Create analysis file path by appending symprec to pred_file basename
    geo_opt_filename = model.geo_opt_path.removesuffix(".json.gz")
    geo_opt_csv_path = f"{geo_opt_filename}-{symprec_str}-{moyo_version}.csv.gz"

    if os.path.isfile(geo_opt_csv_path):
        print(f"{prog_str} {model_label} already analyzed at {geo_opt_csv_path}")
        continue

    # Load model structures
    try:
        if ml_relaxed_structs_path.endswith((".json", ".json.gz", ".json.xz")):
            df_ml_structs = pd.read_json(ml_relaxed_structs_path).set_index(Key.mat_id)
        else:
            raise ValueError(
                "Relaxed structure analysis currently only supports pymatgen JSON, "
                f"got {ml_relaxed_structs_path}"
            )
    except Exception as exc:
        exc.add_note(f"{model_label=} {ml_relaxed_structs_path=}")
        raise

    if debug_mode:
        df_ml_structs = df_ml_structs.head(debug_mode)

    struct_col = geo_opt_metrics.get("struct_col")
    if struct_col not in df_ml_structs:
        struct_cols = [col for col in df_ml_structs if Key.structure in col]
        print(
            f"⚠️ {struct_col=} not found in {model_label}-relaxed structures loaded "
            f"from {ml_relaxed_structs_path}. Did you mean one of {struct_cols}?"
        )
        continue

    # Convert structures
    model_structs = {
        mat_id: Structure.from_dict(struct_dict)
        for mat_id, struct_dict in df_ml_structs[struct_col].items()
    }

    # Analyze symmetry
    model_analysis = analyze_symmetry(
        model_structs,
        pbar=dict(desc=f"{prog_str} Analyzing {model_label} symmetries"),
        symprec=symprec,
    )

    # Compare with DFT reference
    df_ml_geo_analysis = pred_vs_ref_struct_symmetry(
        model_analysis,
        dft_analysis,
        model_structs,
        dft_structs,
        pbar=dict(desc=f"{prog_str} Comparing DFT vs {model_label} symmetries"),
    )

    # Save model results
    df_ml_geo_analysis.to_csv(geo_opt_csv_path)
    print(f"{prog_str} Completed {model_label} and saved results to {geo_opt_csv_path}")

    # Calculate metrics and write to YAML
    df_metrics = geo_opt.calc_geo_opt_metrics(df_ml_geo_analysis)
    geo_opt.write_geo_opt_metrics_to_yaml(df_metrics, model, symprec)

print("All models processed!")
